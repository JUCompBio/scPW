# scGSVA on the Seurat base image, pinned to a commit.
#
# No official scGSVA image exists. Installing from GitHub inside every run would make
# each benchmark depend on the network and on whatever the default branch happens to be
# that day, which is not acceptable in a methods comparison. So the install happens once,
# here, at a fixed commit recorded in the image.
#
# The Seurat base is reused rather than a plain r-base so that Matrix, cluster and the
# rest of the toolchain match the Seurat arm exactly. Any difference between the scGSVA
# and Seurat rows is then attributable to the method, not to a different BLAS or a
# different Matrix version.
#
# The first build attempts failed in ways worth recording, because both errors surfaced
# several packages away from their cause.
#
# Attempt 1:
#   ERROR: dependency 'nloptr' is not available for package 'lme4'
#   ERROR: dependency 'lme4' ... 'car' ... 'rstatix' ... 'scGSVA'
# nloptr builds its bundled NLopt with CMake, which this image does not ship.
#
# Attempt 3, after pinning CRAN, and attempt 4 after rebuilding Matrix from source:
#   lme4CholmodDecomposition.h: error: 'cholmod_factor' does not name a type
#   lme4CholmodDecomposition.h: error: expected template-name before '<' token
# lme4 links against the C-level CHOLMOD headers of Matrix and the Eigen headers of
# RcppEigen, and must be compiled against the same versions it will run with. This image
# carries Matrix 1.6.1.1 and RcppEigen 0.3.3.9.4, both newer than the 2023 lme4 expects,
# and rebuilding Matrix alone does not resolve the Eigen half of the conflict.
#
# The resolution is to stop compiling that chain at all. This image is Ubuntu 20.04
# (focal) with R 4.2.0, and Posit publishes precompiled binaries for exactly that
# combination, built against a mutually consistent set of Matrix and RcppEigen. Binary
# packages cannot exhibit an ABI mismatch against their own build environment, so the
# failure disappears rather than being worked around.
#
# Two details make the binary repository actually serve binaries:
#   * the URL needs the __linux__/focal prefix
#   * HTTPUserAgent must carry the R version, or the server falls back to source tarballs
#     and the compile failures return
#
# Note on why this chain matters at all: rstatix is in scGSVA's Imports for three
# functions (anova_test, t_test, wilcox_test, all in misc.R, none on the scoring path),
# and car, doBy, pbkrtest and lme4 are referenced zero times anywhere in the package. The
# whole chain is transitive, which is why it is worth satisfying cheaply rather than
# fighting.

FROM satijalab/seurat:5.0.0

# Pinned: scGSVA is described by its author as under active development, so an unpinned
# install would not be reproducible.
ARG SCGSVA_COMMIT=cd5a2b1358012b23060b5c6e892c4eb556c65a8b
# Posit's date-stamped binary mirror for Ubuntu 20.04. Date-pinned for reproducibility;
# 2023-04-21 is late in R 4.2's life, so packages are recent but built against this R.
ARG PPM_BINARY=https://packagemanager.posit.co/cran/__linux__/focal/2023-04-21

# Kept for any package that still has to build from source, and for system libraries the
# binaries link against at run time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake libnlopt-dev libgsl-dev libglpk-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Both lines are required. Without the user agent, Posit cannot tell which R version is
# asking and serves source tarballs, which reintroduces the lme4 compile failure.
#
# The timeout matters too. R defaults to 60 seconds per download, and org.Hs.eg.db is an
# 84 MB annotation tarball from bioconductor.org; attempt 5 failed on it after 483 s with
# no compile error, which is the signature of a download timing out rather than anything
# structural. One hour is generous enough that a slow link is not mistaken for a failure.
RUN printf '%s\n' \
    "options(repos = c(CRAN = '${PPM_BINARY}'))" \
    "options(timeout = 3600)" \
    "options(HTTPUserAgent = sprintf('R/%s R (%s)', getRversion(), paste(getRversion(), R.version\$platform, R.version\$arch, R.version\$os)))" \
    >> /usr/local/lib/R/etc/Rprofile.site

# Confirm binaries are actually being served before installing anything heavy: a source
# install of lme4 here is the failure mode this whole layer exists to avoid.
RUN R -q -e "ap <- available.packages(); stopifnot('lme4' %in% rownames(ap)); \
             cat('repo:', getOption('repos')[['CRAN']], '\n'); \
             cat('lme4 offered:', ap['lme4','Version'], '\n')"

RUN R -q -e "install.packages(c('remotes','BiocManager'))"

# Register the Bioconductor repositories alongside the binary CRAN mirror, so dependency
# resolution covers Bioconductor as well. CRAN stays first so its prebuilt binaries win
# over Bioconductor's source copies.
RUN R -q -e "bio <- BiocManager::repositories(); \
             bio <- bio[names(bio) != 'CRAN']; \
             line <- sprintf(\"options(repos = c(CRAN = '%s', %s))\", \
                             '${PPM_BINARY}', \
                             paste(sprintf(\"%s = '%s'\", names(bio), bio), collapse=', ')); \
             cat(line, file='/usr/local/lib/R/etc/Rprofile.site', append=TRUE, sep='\n'); \
             cat(line, '\n')"

# scGSVA declares 23 dependencies and the Seurat image already carries 20 of them. Only
# these three are missing, determined by reading the DESCRIPTION at the pinned commit and
# checking it against installed.packages() rather than by listing packages from memory,
# which is how GO.db was overlooked on the previous attempt.
#
# rstatix is the one that drags in lme4 -> pbkrtest -> doBy -> car. As a prebuilt binary
# that chain is a download instead of the compile that failed three times.
RUN R -q -e "install.packages(c('msigdbr','rstatix')); \
             for (p in c('msigdbr','rstatix','lme4','car')) \
               if (!requireNamespace(p, quietly=TRUE)) stop(paste(p, 'failed to install'))"

# GO.db and GSVA are small. Kept separate from org.Hs.eg.db so that a failure on the
# large annotation download does not discard these, and so the failing step is
# unambiguous.
RUN R -q -e "BiocManager::install(c('GO.db','GSVA'), ask=FALSE, update=FALSE); \
             for (p in c('GO.db','GSVA')) \
               if (!requireNamespace(p, quietly=TRUE)) stop(paste(p, 'failed to install'))"

# org.Hs.eg.db is 84 MB and is the step attempt 5 died on. It is not one of scGSVA's
# declared dependencies; buildAnnot() needs it for human annotation, so it is our
# addition rather than the package's requirement. Its own layer, with the download
# reported, so a repeat failure is clearly a network problem and not a resolution one.
RUN R -q -e "cat('timeout is', getOption('timeout'), 'seconds\n'); \
             BiocManager::install('org.Hs.eg.db', ask=FALSE, update=FALSE); \
             if (!requireNamespace('org.Hs.eg.db', quietly=TRUE)) \
               stop('org.Hs.eg.db failed to install')"

RUN R -q -e "remotes::install_github(paste0('guokai8/scGSVA@', Sys.getenv('SCGSVA_COMMIT')), upgrade='never')"

# Fail the build rather than the run if anything above silently did not take. The first
# attempt produced an image that looked built and had no scGSVA in it.
RUN R -q -e "library(scGSVA); library(GSVA); cat('scGSVA', as.character(packageVersion('scGSVA')), 'OK\n')"

LABEL scgsva_commit=$SCGSVA_COMMIT
LABEL ppm_binary=$PPM_BINARY
