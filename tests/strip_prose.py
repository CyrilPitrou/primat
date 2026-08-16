#!/usr/bin/env python3
"""Prove an edit changed only comments and docstrings, not code.

Writes one comment-stripped copy of every ``primat/``, ``primat-c/`` and
``tests/`` source file into an output directory, mirroring the repo layout.
Snapshot before an editing session and after it: an empty ``diff -r`` between
the two snapshots proves the edit was prose-only, which is what lets a large
comment rewrite skip re-measuring the observables.

Python is parsed to an AST with every docstring node removed and re-dumped via
``ast.dump``, so neither formatting nor comments can survive. C goes through a
comment-stripping state machine that respects string and character literals.
String literals ARE code to this tool: a changed ``CHECK()`` assertion label
shows up as a difference, correctly.

Usage::

    python tests/strip_prose.py . /tmp/before
    ... edit comments ...
    python tests/strip_prose.py . /tmp/after
    diff -rq /tmp/before /tmp/after     # empty => prose-only
"""
import ast
import sys
from pathlib import Path

ROOTS = ("primat", "primat-c", "tests")
_DOCSTRING_OWNERS = (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def strip_py(text: str) -> str:
    """Return an AST dump of ``text`` with every docstring node removed."""
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if not isinstance(node, _DOCSTRING_OWNERS):
            continue
        body = node.body
        if (body and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            # A body of nothing but a docstring still needs a statement.
            node.body = body[1:] or [ast.Pass()]
    ast.fix_missing_locations(tree)
    return ast.dump(tree, indent=1)


def strip_c(text: str) -> str:
    """Return ``text`` with C comments removed and whitespace normalised.

    Quote handling is what makes this trustworthy: ``"/* not a comment */"``
    inside a string literal must survive, and a comment must not be ended by a
    ``*/`` that sits inside one.
    """
    out = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch in ('"', "'"):
            quote = ch
            out.append(ch)
            i += 1
            while i < n:
                if text[i] == "\\":          # escape: consume the pair whole
                    out.append(text[i:i + 2])
                    i += 2
                    continue
                out.append(text[i])
                if text[i] == quote:
                    i += 1
                    break
                i += 1
            continue
        if ch == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                i += 1
            i += 2
            out.append(" ")
            continue
        if ch == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                if text[i] == "\\" and i + 1 < n and text[i + 1] == "\n":
                    i += 2          # a line continuation extends a // comment
                    continue
                i += 1
            out.append(" ")
            continue
        out.append(ch)
        i += 1
    return "\n".join(" ".join(line.split())
                     for line in "".join(out).splitlines() if line.split())


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__.split("Usage::")[1].strip(), file=sys.stderr)
        return 2
    root, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    count = 0
    for r in ROOTS:
        for path in sorted((root / r).rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix == ".py":
                fn = strip_py
            elif path.suffix in (".c", ".h"):
                fn = strip_c
            else:
                continue
            dest = out / path.relative_to(root)
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                dest.write_text(fn(path.read_text()))
            except SyntaxError as exc:
                print(f"SYNTAX ERROR {path.relative_to(root)}: {exc}", file=sys.stderr)
                return 1
            count += 1
    print(f"stripped {count} files -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
