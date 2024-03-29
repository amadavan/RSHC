licenses(["notice"])

# Copies environment to tmp and runs ./configure so we can get the
# configuration header file
genrule(
  name = "configure_hdf5",
  srcs = glob(["**/*"], exclude=["src/H5pubconf.h", "src/libhdf5.settings", "src/libsettings.c"]),
  outs = ["src/H5pubconf.h", "src/libhdf5.settings"],
  cmd = "pushd external/hdf5/; workdir=$$(mktemp -d -t tmp.XXXXXXXXXX); cp -r * $$workdir; pushd $$workdir; ./configure --enable-build-mode=production --disable-hl --with-zlib; popd; popd; cp $$workdir/src/H5pubconf.h $$workdir/src/libhdf5.settings $(@D)/src; rm -rf $$workdir;",
  tools = ["configure"],
  message = "Configuring HDF5",
  visibility = ["//visibility:private"],
)

cc_library(
  name = "h5_includes",
  hdrs = glob(["src/*.h"], exclude=["src/H5pubconf.h"]) + [":configure_hdf5"],
  includes = ["src"],
  visibility = ["//visibility:private"],
)

cc_binary(
  name = "h5make_libsettings",
  srcs = ["src/H5make_libsettings.c"],
  deps = [":h5_includes"],
  visibility = ["//visibility:private"],
)

cc_binary(
  name = "h5detect",
  srcs = ["src/H5detect.c"],
  deps = [":h5_includes"],
  visibility = ["//visibility:private"],
)

genrule(
  name = "create_libsettings",
  srcs = [":configure_hdf5"],
  outs = ["src/libsettings.c"],
  cmd = "cp $(@D)/libhdf5.settings .; $(location //:h5make_libsettings) > $(@D)/libsettings.c; rm libhdf5.settings",
  tools = [":h5make_libsettings"],
  visibility = ["//visibility:private"],
)

genrule(
  name = "create_native",
  outs = ["src/native.c"],
  cmd = "$(location :h5detect) > $(@D)/native.c",
  tools = [":h5detect"],
  visibility = ["//visibility:private"],
)

cc_library(
  name = "hdf5",
  srcs = glob(["src/*.c"], exclude=["src/H5make_libsettings.c", "src/H5detect.c"]) + [":create_native"] + [":create_libsettings"],
  hdrs = glob(["src/*.h"], exclude=["src/H5pubconf.h"]) + [":configure_hdf5"],
  linkopts = ["-ldl"],
  includes = ["src"],
  deps = ["@zlib//:zlib"],
  visibility = ["//visibility:public"],
)

cc_library(
  name = "hdf5_cpp",
  srcs = glob(["c++/src/*.cpp"]),
  hdrs = glob(["c++/src/*.h"]),
  includes = ["c++/src"],
  deps = [":hdf5"],
  visibility = ["//visibility:public"],
)