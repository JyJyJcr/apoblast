{
  stdenvNoCC,
  lib,
  corefonts,
  resvg,
}:

stdenvNoCC.mkDerivation {
  pname = "apoblast-logo";
  version = "1.0.0";

  src = ./.;

  nativeBuildInputs = [
    resvg
  ];

  dontConfigure = true;

  buildPhase = ''
    runHook preBuild

    mkdir -p build

    for svg_file in *.svg; do
      [ -e "$svg_file" ] || continue
      filename=$(basename "$svg_file" .svg)

      cp "$filename.svg" "build/$filename.svg"

      echo "Convert $filename.svg to $filename.png"
      resvg \
        --skip-system-fonts \
        --use-fonts-dir ${corefonts}/share/fonts/truetype \
        -z 4 \
        "build/$filename.svg" \
        "build/$filename.png"
    done

    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall

    mkdir -p "$out/share/apoblast"

    cp build/*.svg "$out/share/apoblast/"
    cp build/*.png "$out/share/apoblast/"

    runHook postInstall
  '';

  meta = {
    description = "SVG and PNG logo for the apoblast project";
    license = lib.licenses.wtfpl;
    platforms = lib.platforms.all;
  };
}
