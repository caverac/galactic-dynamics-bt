window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    macros: {
      // Custom astronomy/physics macros
      rmarcsec: "\\mathrm{arcsec}",
      rmAU: "\\mathrm{AU}",
      rmcentury: "\\mathrm{century}",
      rmpc: "\\mathrm{pc}",
      rmkpc: "\\mathrm{kpc}",
      rmMpc: "\\mathrm{Mpc}",
      rmyr: "\\mathrm{yr}",
      rmMyr: "\\mathrm{Myr}",
      rmGyr: "\\mathrm{Gyr}",
      rmkm: "\\mathrm{km}",
      rmms: "\\mathrm{m\\,s}",
      rmkms: "\\mathrm{km\\,s}",
      // Solar system units
      Msun: "{\\mathrm{M}}_{\\odot}",
      Lsun: "L_{\\odot}",
      Rsun: "R_{\\odot}",
    },
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

document$.subscribe(() => {
  MathJax.typesetPromise();
});
