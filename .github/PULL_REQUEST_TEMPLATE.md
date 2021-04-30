## Overview

Short description of what the PR addresses, eg. new feature or bug fix.

## PR checklist
 - [ ] Adds documentation for the new scripts, functions you added
 - [ ] If scripts are inccluded in the docker container, make sure to build and push to DockerHub
 - [ ] Link issues that this PR addresses from the `Linked issues` section in the bottom right of the Pull Request page (after the PR is created)
 - [ ] Rebuild and push the container if you have updated the Dockerfile
 - [ ] Successfully completed job pointing to latest commit in PR
 - [ ] Example data added in S3 in the following location s3://lifebit-featured-datasets/pipelines/<name-of-repo>/testdata/ 
> NOTE: (ask [@cgpu](https://github.com/cgpu) if you don't have access to the S3 bucket)