"use strict";
const form = document.getElementById("login");
const err = document.getElementById("err");
const go = document.getElementById("go");
fetch("/api/session").then(r => r.json()).then(s => {
  if (s.user) location.replace("/");
  else if (!s.account) err.textContent = "No password is set yet. On main run: .venv/bin/python -m dashboard passwd";
}).catch(() => {});
form.addEventListener("submit", async (e) => {
  e.preventDefault();
  err.textContent = "";
  go.disabled = true;
  try {
    const r = await fetch("/api/login", {
      method: "POST", headers: {"Content-Type": "application/json", "X-Minos": "1"},
      body: JSON.stringify({user: form.user.value, password: form.password.value}),
    });
    if (r.ok) { location.replace("/"); return; }
    const body = await r.json().catch(() => ({}));
    err.textContent = body.error || `Login failed (${r.status})`;
    form.password.value = "";
    form.password.focus();
  } catch (x) {
    err.textContent = "The server did not answer. Is the dashboard service running?";
  } finally {
    go.disabled = false;
  }
});
