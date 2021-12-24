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
    rsync
  ] ++ lib.optional enableDoxygen [ doxygen graphviz ];

  buildInputs = with pkgs; [
    openblas
    lapack-reference
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
    mkdir -p $out/bin
    cp itpp-config $out/bin
    chmod a+rx $out/bin/itpp-config

    mkdir -p $out/lib
    cp itpp/libitpp*.dylib $out/lib

    # Headers
    mkdir -p $out/include/itpp
    rsync --recursive --prune-empty-dirs --exclude='config_msvc.h' --include='*.h' --include='*/' --exclude='*' $src/itpp/ $out/include/itpp/
    cp itpp/config.h $out/include/itpp
    cp itpp/itexports.h $out/include/itpp

    # Extra (Matlab, Python)
    mkdir -p $out/share/itpp
    cp $src/extras/itsave.m $out/share/itpp/
    cp $src/extras/itload.m $out/share/itpp/
    cp $src/extras/pyitpp.py $out/share/itpp/
    cp $src/extras/gdb_macros_for_itpp $out/share/itpp/

    # Pkg-config support
    mkdir -p $out/lib/pkgconfig
    cp itpp.pc $out/lib/pkgconfig

    # Man page
    mkdir -p $out/share/man/man1
    cp itpp-config.1 $out/share/man/man1

    # Doxygen documentation
    if [ -d "html" ]; then
      mkdir -p $out/share/doc/itpp
      cp -r html $out/share/doc/itpp
      cp $src/doc/images/itpp_logo.png $out/share/doc/itpp/html
    fi
  '';
}
