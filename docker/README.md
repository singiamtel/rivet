# DockerSetup.md

This is a guide to building Rivet Docker images. Docker is used to
create consistent, isolated environments that can be easily shared and
deployed. This covers the steps necessary to configure your system for
Docker, build and push Docker images, and integrate Docker into your
CI/CD pipeline.


## Prerequisites

Before you begin, ensure that your system meets the following requirements:

- **Docker**: Installed on your system. Follow the [official Docker installation guide](https://docs.docker.com/engine/installation/) if needed.
- **QEMU**: Required for multi-architecture builds.
- **Git**: To clone the Rivet repository and manage source code.



## Installation and Setup

First, obviously, clone the Rivet git repo.


### Multi-architecture builds

To set up your environment for multi-architecture builds, run the following:

```bash
chmod +x ./docker/setup-multiarch-builder.sh
./docker/setup-multiarch-builder.sh
```

This script installs QEMU, sets up Docker Buildx, and configures your
system for building images for multiple architectures.


### Common Docker setup

The `docker-common.sh` script contains common setup logic used across
different Docker build scripts. It is automatically sourced by
`build.sh` and handles tasks like setting up build flags and the
Docker build function. You don't need to run this, but any updates to
general configs that affect all image builds should go here.


## Building Docker images

To build Docker images, use the `build.sh` script in each image directory:

```bash
cd rivet/
FLAG1=foo FLAG2=bar ./build.sh
```

The flags given above are just environment variables passed to the script.
The key ones are:

- **MULTIARCH=1**: Enables multi-architecture builds.
- **PUSH=1**: Pushes the built images to Docker Hub.
- **FORCE=1**: Forces a rebuild of the images, ignoring any cached layers.
- **SLEEP=30s**: Leave a gap of the specified time between image pushes, to avoid overload/rejection

### Push to Docker Hub

If you want to push the images to Docker Hub, ensure that the `PUSH`
flag is set to `1`:

```bash
PUSH=1 ./build.sh
```

### Tagging and Versioning

The `build.sh` script handles tagging during the build process. Tags
are applied based on the branch name and architecture, ensuring that
your images are correctly versioned and identifiable.


## CI/CD integration

The `.gitlab-ci.yml` file is configured to be able to build and test images.

### CI Configuration

The CI/CD pipeline is set up to run the `build.sh` script with the
appropriate flags for multi-architecture builds and pushing images to
Docker Hub. Here’s a breakdown of the relevant sections:

```yaml
stages:
  - build

variables:
  DOCKER_DRIVER: overlay2
  IMAGE_NAME: rivet
  FORCE: "1"
  PUSH: "1"
  MULTIARCH: "1"

build_docker_image:
  stage: build
  image: docker:latest
  services:
    - docker:dind
  before_script:
    - chmod +x ./docker/setup-multiarch-builder.sh
    - ./docker/setup-multiarch-builder.sh
  script:
    - echo "$DOCKER_HUB_PASSWORD" | docker login -u "$DOCKER_HUB_USERNAME" --password-stdin
    - chmod +x ./docker/rivet/build.sh
    - cd ./docker/rivet
    - ./build.sh
  only:
    - feature/multiarch-docker-images-from-ci
```

### Running the pipeline

To trigger the CI/CD pipeline, push your changes to the repository or
create a merge request targeting the
`feature/multiarch-docker-images-from-ci` branch.

TODO: change this to trigger on release commits or manual build scheduling.
