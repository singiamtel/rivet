#! /usr/bin/env bash

## Exit immediately on command errors
set -e

## Determine platform flags for multi-architecture builds
if [[ "$MULTIARCH" = 1 ]]; then
    PLATFLAGS="--platform linux/amd64,linux/arm64"
else
    PLATFLAGS=""
fi

BUILDFLAGS="$PLATFLAGS $DOCKERFLAGS"

## Function to wrap Docker commands
function xdocker { echo "docker $@"; docker "$@"; }

## Define the build function based on the MULTIARCH flag
if [[ "$MULTIARCH" = 1 ]]; then
    if [[ "$PUSH" = 1 ]]; then PUSH="--push"; fi
    function dx_build { xdocker buildx build -f ./Dockerfile --progress=plain "$@" $PUSH .; }
else
    function dx_build { xdocker build -f ./Dockerfile --progress=plain "$@" .; }
fi

# ## If Docker Hub username is not found and PUSH is set, prompt the user to enter it
# ## Automatically retrieve the Docker Hub username from the Docker config file
# DOCKER_HUB_USERNAME=${DOCKER_HUB_USERNAME:-$(jq -r '.auths["https://index.docker.io/v1/"].auth' ~/.docker/config.json | base64 --decode | cut -d: -f1)}
# if [[ -z "$DOCKER_HUB_USERNAME" && "$PUSH" = 1 ]]; then
#     read -p "Docker Hub username not found. Please enter your Docker Hub username: " DOCKER_HUB_USERNAME
# fi

## If the FORCE flag is passed, ignore the Docker build cache
if [[ "$FORCE" = 1 ]]; then
    BUILDFLAGS="$BUILDFLAGS --no-cache"
fi
