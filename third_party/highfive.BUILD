licenses(["notice"])

cc_library(
    name = "h5",
    hdrs = glob(["include/highfive/**/*.hpp"]),
    strip_include_prefix = "include",
    includes = ["include"],
    deps = ["@hdf5//:hdf5"],
    visibility = ["//visibility:public"],
)
