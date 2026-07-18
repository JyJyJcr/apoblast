{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
    nix-jyjyjcr = {
      url = "github:jyjyjcr/nix-jyjyjcr";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs =
    {
      nixpkgs,
      flake-utils,
      nix-jyjyjcr,
      ...
    }:
    (flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfreePredicate =
            pkg:
            builtins.elem (nixpkgs.lib.getName pkg) [
              "corefonts"
            ];
        };
        apoblast-logo = pkgs.callPackage ./assets/logo.nix { };

        pkgs-dev = import nixpkgs {
          inherit system;
          overlays = [ nix-jyjyjcr.overlays.default ];
        };
      in
      {
        packages.apoblast-logo = apoblast-logo;

        devShells = pkgs-dev.alt-shell.mkCommonShells { } {
          packages = [
            pkgs-dev.uv
            pkgs-dev.python313
          ];
        };
      }
    ))
    // {
      overlays.default = final: prev: {
        apoblast-logo = prev.callPackage ./assets/logo.nix { };
      };
    };
}
