workspace(name = "rshc")

# Local configuration
load("@rshc//third_party/systemlibs:syslibs_configure.bzl", "syslibs_configure")
load("@rshc//third_party/py:python_configure.bzl", "python_configure")
load("@rshc//third_party:repo.bzl", "rshc_http_archive")

syslibs_configure(name = "local_config_syslibs")
python_configure(name = "local_config_python")

# Local libraries
load("@rshc//third_party/gurobi:build_defs.bzl", "gurobi_repository")

gurobi_repository(
        name = "gurobi",
        build_file = "@rshc//third_party/gurobi:gurobi.BUILD",
)

# External libraries
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

# Rule repository
http_archive(
   name = "rules_foreign_cc",
   strip_prefix = "rules_foreign_cc-0.9.0",
   sha256 = "2a4d07cd64b0719b39a7c12218a3e507672b82a97b98c6a89d38565894cf7c51",
   urls = [
       "https://github.com/bazelbuild/rules_foreign_cc/archive/0.9.0.tar.gz",
       ]
)

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")
rules_foreign_cc_dependencies()

# Libraries
http_archive(
    name = "eigen",
    build_file = "@rshc//third_party:eigen.BUILD",
    sha256 = "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72",
    strip_prefix = "eigen-3.4.0",
    urls = [
        "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz",
    ],
)

http_archive(
    name = "phasor",
    build_file = "@rshc//third_party:phasor.BUILD",
    sha256 = "4179d34f9ed7a1add7f1ee348329e664fa41696dca79415e0ea67e365620937a",
    strip_prefix = "Phasor-0.0.1",
    urls = [
        "https://github.com/amadavan/Phasor/archive/refs/tags/v0.0.1.tar.gz"
    ]
)

http_archive(
    name = "phasor-nofs",
    build_file = "@rshc//third_party:phasor.BUILD",
    sha256 = "a049ae5f9ee1dda04a4cd7b45400f46fa27ef8817db1270481d4bbe664892dfb",
    strip_prefix = "Phasor-0.0.1b",
    urls = [
        "https://github.com/amadavan/Phasor/archive/refs/tags/v0.0.1b.tar.gz"
    ]
)

git_repository(
    name = "cxxopts",
    commit = "c74846a891b3cc3bfa992d588b1295f528d43039",
    shallow_since = "1634764013 +1100",
    remote = "https://github.com/jarro2783/cxxopts.git"
)

rshc_http_archive(
    name = "hdf5",
    build_file = "@rshc//third_party:hdf5.BUILD",
    system_build_file = "@rshc//third_party/systemlibs:hdf5.BUILD",
    sha256 = "1a88bbe36213a2cea0c8397201a459643e7155c9dc91e062675b3fb07ee38afe",
    strip_prefix = "hdf5-1.12.2",
    urls = [
        "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.bz2",
        ]
)
rshc_http_archive(
    name = "zlib",
    build_file = "@rshc//third_party:zlib.BUILD",
    system_build_file = "@rshc//third_party/systemlibs:zlib.BUILD",
    sha256 = "b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30",
    strip_prefix = "zlib-1.2.13",
    urls = [
        "https://github.com/madler/zlib/releases/download/v1.2.13/zlib-1.2.13.tar.gz",
        "https://zlib.net/zlib-1.2.13.tar.gz",
    ],
)

rshc_http_archive(
    name = "szip",
    build_file = "@rshc//third_party:szip.BUILD",
    system_build_file = "@rshc//third_party/systemlibs:szip.BUILD",
    sha256 = "21ee958b4f2d4be2c9cabfa5e1a94877043609ce86fde5f286f105f7ff84d412",
    strip_prefix = "szip-2.1.1",
    urls = [
        "https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz",
    ]
)

http_archive(
    name = "highfive",
    build_file = "@rshc//third_party:highfive.BUILD",
    sha256 = "6826471ef5c645ebf947d29574b302991525a8a8ff1ef687aba7311d9a0ea36f",
    strip_prefix = "HighFive-2.4.1",
    urls = [
        "https://github.com/BlueBrain/HighFive/archive/refs/tags/v2.4.1.tar.gz"
    ]
)

http_archive(
    name = "mosek_linux",
    build_file = "@rshc//third_party/mosek:mosek.BUILD",
    sha256 = "edfd81be24a6d1dca34b97a783d009167467092057afca4e3555d0b97305bbcc",
    strip_prefix = "mosek/9.3/tools/platform/linux64x86",
    urls = [
        "https://download.mosek.com/stable/9.3.14/mosektoolslinux64x86.tar.bz2"
    ]
)

http_archive(
    name = "mosek_darwin",
    build_file = "@rshc//third_party/mosek:mosek.BUILD",
    sha256 = "f843aeed445c3f30061120c8837c5187dd7cf0b6307712ad99fb0ee480f71f4b",
    strip_prefix = "mosek/9.3/tools/platform/osx64x86",
    urls = [
        "https://download.mosek.com/stable/9.3.14/mosektoolsosx64x86.tar.bz2"
    ]
)

http_archive(
    name = "mosek_windows",
    build_file = "@rshc//third_party/mosek:mosek.BUILD",
    sha256 = "e2fc2b460f632fe3806bdc74d4329336d8284bfe7043190c18ae998f7caaa2aa",
    strip_prefix = "mosek/9.3/tools/platform/win64x86",
    urls = [
        "https://download.mosek.com/stable/9.3.14/mosektoolswin64x86.zip"
    ]
)

# Compilation database (for vscode intellisense)
http_archive(
    name = "com_grail_bazel_compdb",
    strip_prefix = "bazel-compilation-database-0.5.2",
    sha256 = "d32835b26dd35aad8fd0ba0d712265df6565a3ad860d39e4c01ad41059ea7eda",
    urls = ["https://github.com/grailbio/bazel-compilation-database/archive/0.5.2.tar.gz"],
)

load("@com_grail_bazel_compdb//:deps.bzl", "bazel_compdb_deps")
bazel_compdb_deps()