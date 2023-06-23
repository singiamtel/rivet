#! /usr/bin/env bash
# Wrap the bash setup so it'll still work as a Singularity container

. /etc/profile
. /etc/bash.bashrc
. /herwig/bin/activate
for i in /etc/profile.d/*; do source $i; done

if [[ -n "$@" ]]; then
    eval "$@"
else
    bash --rcfile /etc/bash.bashrc
fi
