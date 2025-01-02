test "$LATEST" = 1 && TAGS+=" latest"
for tag in $TAGS; do TAGFLAGS+=" -t $PKG:$tag"; done

dx_build $BUILDFLAGS $TAGFLAGS

if [[ "$PUSH" = 1 ]]; then
    for tag in $TAGS; do xdocker push $PKG:$tag; done
fi
