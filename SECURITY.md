# Security policy

## Supported versions

Only the latest released version receives fixes.

## Reporting a vulnerability

Please do **not** open a public issue. Use GitHub's
[private vulnerability reporting](../../security/advisories/new) or write to **efirvida@gmail.com**
with the affected version, a description and, if possible, a minimal reproduction. You will get an
answer within a few days.

## Scope notes

* The geometry exporters render Jinja2 templates and write files to a directory chosen by the caller;
  treat model names and output paths from untrusted sources accordingly.
* `gearbox.fem` calls the external programs Gmsh and CalculiX; their own advisories apply.
* **Engineering results are not a security matter.** A wrong stress or safety factor is a bug: report it
  as an [issue](../../issues/new/choose) with the standard clause. The results do not replace a design
  review by a qualified engineer.
