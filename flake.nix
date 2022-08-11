{
  description = "VAMPyR";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    pypi-deps-db = {
      url = "github:DavHau/pypi-deps-db";
      flake = false;
    };
    mach-nix = {
      url = "mach-nix/3.4.0";
      inputs.pypi-deps-db.follows = "pypi-deps-db";
    };
  };

  outputs = { self, nixpkgs, flake-utils, pypi-deps-db, mach-nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        pythonEnv = mach-nix.lib."${system}".mkPython {
          requirements = builtins.readFile ./requirements.txt + ''
                  # additional dependencies for local work
                  #jupyterlab
            pre-commit
          '';
        };
      in
      {
        devShell = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            cmake
            eigen
            gcc
            openmpi
            ninja
          ];
          buildInputs = [
            pythonEnv
          ];
          NINJA_STATUS = "[Built edge %f of %t in %e sec] ";
        };
      });
}
