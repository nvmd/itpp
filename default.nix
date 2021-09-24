{ pkgs ? import <nixpkgs> {}
, enableDoxygen ? false
}:

let gtestSrc = pkgs.fetchFromGitHub {
      owner = "google";
      repo = "googletest";
      rev = "release-1.11.0";
      hash = "sha256-SjlJxushfry13RGA7BCjYC9oZqV4z6x8dOiHfl/wpF0=";
    };

in pkgs.stdenv.mkDerivation rec {
  name = "itpp";
  # version = "4.3.1";

  src = builtins.path { name = "itpp"; path = ./.; };

  enableParallelBuilding = true;

  nativeBuildInputs = with pkgs; [
    cmake
    python
  ] ++ lib.optional enableDoxygen [ doxygen graphviz ];

  buildInputs = with pkgs; [
    blas
    lapack
    fftw
    gtestSrc
  ];

  configurePhase = ''
    mkdir -p build && cd build
    cmake -DOLD_TESTS=on -DGTEST_DIR=${gtestSrc}/googletest ..
  '';

  buildPhase = ''
    make
  '';

  doCheck = true;
  checkPhase = ''
    echo Running old unit tests
    python ../extras/check_tests.py -r ../tests -w tests
    echo Running GTest-based unit tests
    ./gtests/itpp_gtests
  '';

  installPhase = ''
    mkdir -p $out $out/build
    cp -r $src/* $out/
    cp -r . $out/build/
    echo installPhase is not yet implemented
  '';
}
