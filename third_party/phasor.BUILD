
PHASOR_HEADER_FILES = glob(["include/**"])

cc_library(
    name = "phasor",
    hdrs = PHASOR_HEADER_FILES,
    includes = ["include"],
    deps = [
        "@eigen//:eigen",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "phasor_header_files",
    srcs = PHASOR_HEADER_FILES,
    visibility = ["//visibility:public"]
)