{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell rec {

   nativeBuildInputs = with pkgs; [
    cmake
    python
    graphviz
    doxygen
    rsync

    # debugging tools
    gdb
   ];

   buildInputs = with pkgs; [
    blas
    lapack
    fftw
   ];
}