#!/bin/sh

set -e

apk add --no-cache curl bash qemu qemu-system-arm qemu-system-aarch64 qemu-system-riscv64

mkdir -p ~/.docker/cli-plugins

curl -fsSL https://github.com/docker/buildx/releases/download/v0.7.1/buildx-v0.7.1.linux-amd64 -o ~/.docker/cli-plugins/docker-buildx
chmod +x ~/.docker/cli-plugins/docker-buildx

docker buildx create --use --name multiarch_builder

# Optional: Bootstrap the Buildx builder (required for some environments)
docker buildx inspect --bootstrap

