"use strict";
// Runs before the page draws (loaded without defer): the theme saved by the toggle, so there is no flash.
(function () {
  try {
    const t = localStorage.getItem("minos-theme");
    if (t === "dark" || t === "light") document.documentElement.dataset.theme = t;
  } catch (e) { /* storage blocked: follow the system setting */ }
})();
