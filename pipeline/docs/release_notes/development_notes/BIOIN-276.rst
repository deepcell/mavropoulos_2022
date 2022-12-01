=====================================
Create Bitbucket pipeline for askCell
=====================================

Description:

This update creates bitbucket CI / CD pipeline for the repo. The pipeline:

 * is triggered on pull-request

 * uses `Dockerfile` that is already created in the repo, rather than create a duplicate of what is written in `Dockerfile` in `bitbucket-pipelines.yml`

 * caches docker images

Stories:

    1. `BIOIN-276 <https://deepcellbio.atlassian.net/browse/BIOIN-276>`_: Create Bitbucket pipeline for askCell
