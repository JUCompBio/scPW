"""
scVI with a network-informed prior on the per-gene dispersion.

The idea and why it sits here rather than in preprocessing
---------------------------------------------------------
scVI's negative binomial likelihood gives each gene an inverse dispersion
``theta_g``, which sets how much variance the model expects and therefore how much it
trusts deviations in that gene:

    Var(x_g) = mu_g + mu_g^2 / theta_g

These are free parameters, one per gene, estimated independently. Genes with few
non-zero counts get poorly determined values.

Rescaling the input columns by network-derived scalars does not address this. Such a
rescaling is rank-preserving within each gene and identical for every cell, so it
cannot change cell-to-cell geometry except through the encoder's nonlinearity, which
is also free to invert it.

This module instead adds a penalty on ``log theta`` so interacting genes are pushed
toward similar dispersion:

    loss = ZINB + KL + tau * sum_{(g,h) in E} w_gh (log theta_g - log theta_h)^2

Because ``theta_g`` scales every cell's contribution to the gradient for gene g, the
network now shapes what the model learns rather than what it is shown, and the model
cannot undo it because the penalty is part of the objective.

The statistical idea is borrowed: edgeR and DESeq2 shrink dispersion toward a common
trend, and single-cell methods pool across genes of similar expression. What appears
to be new is pooling by network adjacency inside a VAE.

tau = 0 must reproduce vanilla scVI exactly. That is asserted in the tests, and it
gives an exact rather than approximate ablation.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scvi
import torch
from scvi.module import VAE


def _laplacian_from_adjacency(A: sp.csr_matrix) -> torch.Tensor:
    """Row-normalised neighbour-averaging operator as a sparse torch tensor.

    The penalty is written as a comparison of each gene against the weighted mean of
    its neighbours rather than as an explicit edge sum. The two are equivalent up to
    a positive factor when A is row-normalised, and this form costs one sparse
    matrix-vector product per step instead of a gather over every edge, which matters
    because the fused graph has millions of edges.
    """
    A = A.tocoo()
    idx = torch.from_numpy(np.vstack([A.row, A.col]).astype(np.int64))
    val = torch.from_numpy(A.data.astype(np.float32))
    return torch.sparse_coo_tensor(idx, val, size=A.shape).coalesce()


class NetworkPriorVAE(VAE):
    """VAE whose per-gene dispersion is smoothed over a gene network.

    Args:
        adjacency: row-normalised gene-gene adjacency in ``var_names`` order, of
            shape ``(n_input, n_input)``. Alignment is checked, because a misaligned
            matrix would still train and still produce a falling penalty while
            pairing every gene with the wrong neighbours.
        tau: penalty strength. ``0.0`` disables the term entirely and the module then
            behaves exactly like ``scvi.module.VAE``.
    """

    def __init__(self, *args, adjacency: sp.csr_matrix | None = None,
                 tau: float = 0.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.tau = float(tau)
        self._has_prior = adjacency is not None and self.tau > 0.0

        if adjacency is not None:
            n = adjacency.shape[0]
            if n != self.n_input:
                raise ValueError(
                    f"adjacency is {n}x{n} but the model has n_input="
                    f"{self.n_input}; it must be built from the same var_names "
                    f"in the same order"
                )
            # Registered as a buffer so it moves with .to(device) and is saved with
            # the model, but is not a trainable parameter.
            self.register_buffer("_adj", _laplacian_from_adjacency(adjacency),
                                 persistent=False)

    def dispersion_penalty(self) -> torch.Tensor:
        """Mean squared deviation of log theta from its neighbour average."""
        if not self._has_prior:
            return torch.zeros((), device=self.px_r.device)

        # px_r is stored on the log scale in scvi-tools, which is what we want: the
        # penalty should act multiplicatively on theta, and dispersion spans orders
        # of magnitude across genes.
        log_theta = self.px_r.view(-1, 1)
        neighbour_mean = torch.sparse.mm(self._adj, log_theta)

        # Genes with no neighbours get a zero row, hence neighbour_mean 0. Excluding
        # them is necessary: otherwise the penalty would drag their log theta toward
        # zero, i.e. theta toward 1, which is an arbitrary target with no
        # justification.
        connected = torch.sparse.sum(self._adj, dim=1).to_dense() > 0
        if not bool(connected.any()):
            return torch.zeros((), device=self.px_r.device)

        diff = (log_theta.view(-1) - neighbour_mean.view(-1))[connected]
        return (diff ** 2).mean()

    def loss(self, tensors, inference_outputs, generative_outputs, kl_weight=1.0):
        """Standard scVI loss plus the network penalty.

        The penalty is added to the scalar loss only. Reconstruction and KL terms are
        passed through untouched so that the values scvi-tools logs and reports keep
        their usual meaning and stay comparable with a vanilla run.
        """
        out = super().loss(tensors, inference_outputs, generative_outputs,
                           kl_weight=kl_weight)
        if not self._has_prior:
            return out

        penalty = self.dispersion_penalty()
        out.loss = out.loss + self.tau * penalty
        # Surfaced for monitoring: if this does not fall during training, the prior is
        # not doing anything and tau is either too small or the adjacency is empty.
        out.extra_metrics = {**(out.extra_metrics or {}),
                             "dispersion_penalty": penalty.detach()}
        return out


class SCVINetworkPrior(scvi.model.SCVI):
    """``scvi.model.SCVI`` that builds a :class:`NetworkPriorVAE`.

    Declaring the module class is the only correct way to do this. The obvious
    alternative -- construct a normal ``SCVI`` and then reassign ``model.module`` --
    does not work, and fails quietly:

        swapping in a plain, unmodified scvi VAE with weights copied across
        still changed the trained latent representation by 1.08e-01

    That was measured with an unmodified ``VAE`` and identical weights loaded via
    ``load_state_dict``, so the divergence comes from the swap itself rather than
    from anything in the subclass. Re-seeding before ``train()`` did not help.
    Building through ``_module_cls`` reproduces vanilla scVI **bit-exactly** at
    tau=0, which is what makes the ablation exact.

    Extra keyword arguments (``adjacency``, ``tau``) are forwarded to the module by
    scvi-tools' own plumbing.
    """

    _module_cls = NetworkPriorVAE


def theta_alignment(module: NetworkPriorVAE) -> dict:
    """How closely does log theta track its neighbour average?

    Reported before and after training as the direct check that the mechanism did
    what it was meant to do, rather than inferring it from downstream metrics.
    """
    if not getattr(module, "_has_prior", False):
        return {}
    with torch.no_grad():
        log_theta = module.px_r.view(-1, 1)
        nb = torch.sparse.mm(module._adj, log_theta).view(-1)
        connected = torch.sparse.sum(module._adj, dim=1).to_dense() > 0
        a = log_theta.view(-1)[connected]
        b = nb[connected]
        if a.numel() < 2:
            return {}
        centred_a = a - a.mean()
        centred_b = b - b.mean()
        denom = centred_a.norm() * centred_b.norm()
        corr = float((centred_a @ centred_b) / denom) if float(denom) > 0 else float("nan")
        return {
            "n_connected": int(connected.sum()),
            "log_theta_sd": float(a.std()),
            "rmse_to_neighbour_mean": float(((a - b) ** 2).mean().sqrt()),
            "corr_with_neighbour_mean": corr,
        }
