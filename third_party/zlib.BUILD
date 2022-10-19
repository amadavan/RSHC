licenses(["notice"])

genrule(
    name = "configure_zlib",
    srcs = glob(["**/*"], exclude=["zconf.h"]),
    outs = ["zconf.h"],
    cmd = "pushd external/zlib/; workdir=$$(mktemp -d -t tmp.XXXXXXXXXX); cp -r * $$workdir; pushd $$workdir; ./configure; popd; popd; cp $$workdir/zconf.h $(@D)/zconf.h; rm -rf $$workdir;",
    tools = ["configure"],
    message = "Configuring ZLIB",
    visibility = ["//visibility:private"],
)

cc_library(
    name = "zlib",
    srcs = glob(["*.c"]),
    hdrs = glob(["*.h"]) + [":configure_zlib"],
    local_defines = select({
        "@rshc//:linux_x86_64": [],
        "@rshc//:darwin": ["HAVE_UNISTD_H"],
        "@rshc//:windows": [],
        "//conditions:default": [],
    }),
    includes = ["."],
    visibility = ["//visibility:public"],
)