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

  outputs = { nixpkgs, flake-utils, nix-jyjyjcr, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ nix-jyjyjcr.overlays.default ];
        };
      in {
        devShells = pkgs.alt-shell.mkCommonShells { } {
          packages = [ pkgs.uv pkgs.python313 ];
        };
      });
}
