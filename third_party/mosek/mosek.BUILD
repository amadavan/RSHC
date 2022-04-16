
cc_library(
  name = "mosek_headers",
  hdrs = glob(["h/*.h"]),
  strip_include_prefix = "h",
  visibility = ["//visibility:public"],
)

cc_library(
  name = "mosek_lib_linux",
  srcs = [
    "bin/libmosek64.so.9.3"
  ],
  visibility = ["//visibility:public"],
)

cc_library(
  name = "mosek_lib_windows",
  srcs = [
    "bin/mosek64_9_3.lib"
  ]
)

cc_library(
  name = "mosek_lib_darwin",
  srcs = [
    "bin/libmosek64.dylib"
  ],
  visibility = ["//visibility:public"],
)

cc_library(
  name = "fusion_cxx",
  srcs = glob(["src/fusion_cxx/*.cc"]),
  includes = ["src/fusion_cxx"],
  hdrs = glob(["src/fusion_cxx/*.h"]),
  strip_include_prefix = "src",
  deps = select({
    "@rshc//:linux_x86_64": [
      "@mosek_linux//:mosek_headers",
      "@mosek_linux//:mosek_lib_linux",
    ],
    "@rshc//:windows": [
      "@mosek_windows//:mosek_headers",
      "@mosek_windows//:mosek_lib_windows",
    ],
    "@rshc//:darwin": [
      "@mosek_darwin//:mosek_headers",
      "@mosek_darwin//:mosek_lib_darwin",
    ],
    "//conditions:default": [],
  }),
  visibility = ["//visibility:public"],
)