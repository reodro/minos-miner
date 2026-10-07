"use strict";
const form = document.getElementById("login");
const err = document.getElementById("err");
const go = document.getElementById("go");

fetch("/api/session").then(r => r.json()).then(s => {
  if (s.user) location.replace("/");
  else if (!s.account) err.textContent = "No password is set yet. On main run: .venv/bin/python -m dashboard passwd";
}).catch(() => {});

async function login() {
  err.textContent = "";
  const user = document.getElementById("user").value;
  const password = document.getElementById("password").value;
  if (!password) { err.textContent = "Enter the password."; return; }
  go.disabled = true;
  try {
    const r = await fetch("/api/login", {
      method: "POST", headers: {"Content-Type": "application/json", "X-Minos": "1"},
      body: JSON.stringify({user, password}),
    });
    if (r.ok) { location.replace("/"); return; }
    const body = await r.json().catch(() => ({}));
    err.textContent = body.error || `Login failed (${r.status})`;
    document.getElementById("password").value = "";
    document.getElementById("password").focus();
  } catch (x) {
    err.textContent = "The server did not answer. Is the dashboard service running?";
  } finally {
    go.disabled = false;
  }
}

go.addEventListener("click", login);
form.addEventListener("submit", (e) => { e.preventDefault(); login(); });
form.addEventListener("keydown", (e) => {
  if (e.key === "Enter") { e.preventDefault(); login(); }
});
