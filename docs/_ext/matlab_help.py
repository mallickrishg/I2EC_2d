"""Convert MATLAB help text to reStructuredText for sphinxcontrib-matlabdomain.

Connected to ``autodoc-process-docstring``. The ``.m`` files are never
modified: each help block is parsed as written and rewritten as reST.

Expected layout of a help block (after the leading ``%`` is stripped)::

    NAME  One-line summary that may wrap

      [a, b] = NAME(x, y)            (usage line, dropped)

    Prose paragraphs. Indented lines are equations or code.

    INPUTS
      x      [n x 1] description, continuing
             on indented lines
      opts   struct:
        .tol     tolerance            [1e-6]

    OUTPUTS
      a      [n x m] description

    NOTES / EXAMPLE / SEE ALSO / REFERENCES / AUTHORS / other ALL-CAPS header

Inputs become ``:param:``/``:type:`` fields, outputs a single ``:returns:``
bullet list, other sections rubrics. If a block cannot be parsed, a warning
is logged and the block is shown verbatim.

Set ``MATLAB_HELP_DUMP=<file.json>`` to write every raw block and its reST
to a JSON file (used by ``check_help.py``).
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

try:
    from sphinx.util import logging as _sphinx_logging

    logger = _sphinx_logging.getLogger(__name__)
except ImportError:  # converter used without Sphinx
    import logging

    logger = logging.getLogger(__name__)

# ==========================================================================
# Section headers
# ==========================================================================

_KEYWORDS = [  # (keyword regex, section kind)
    (r"inputs?|arguments\s*:?\s*\(\s*input\s*\)|parameters?", "params"),
    (r"options|name-value(?:\s+(?:pairs?|arguments))?", "params"),
    (r"outputs?|arguments\s*:?\s*\(\s*output\s*\)|returns?", "returns"),
    (r"description", "prose"),
    (r"examples?", "example"),
    (r"see\s+also", "seealso"),
    (r"references?", "refs"),
    (r"author\(s\)|authors?", "author"),
    (r"notes?", "rubric"),
    (r"syntax|usage", "usage"),
    (r"properties", "properties"),
]
_KEYWORD_RE = re.compile(
    r"^(?P<kw>" + "|".join(f"(?:{k})" for k, _ in _KEYWORDS) + r")"
    r"(?P<num>\s*#\s*\d+)?(?P<sep>\s*[:\-](?=\s|$))?(?P<rest>.*)$",
    re.IGNORECASE,
)
_CAPS_WORD = r"(?:[A-Z][A-Z0-9_()/\-]*[A-Z0-9)]|[A-Z]|&)"
_GENERIC_HEADER_RE = re.compile(
    rf"^(?P<title>{_CAPS_WORD}(?:[ \t]+{_CAPS_WORD})*)(?P<colon>\s*:)?(?P<rest>.*)$"
)
_BANNER_RE = re.compile(r"^\s*(?:[%=\-*#~_]\s*){4,}$")
_IDENT_LIST_RE = re.compile(r"[\w.]+(?:\s*,?\s*[\w.]+)*\s*[.,]?")


def _classify_keyword(kw: str) -> str:
    for pattern, kind in _KEYWORDS:
        if re.fullmatch(pattern, kw, re.IGNORECASE):
            return kind
    return "rubric"


def _title(text: str) -> str:
    text = re.sub(r"\s+", " ", text.strip().rstrip(":").strip())
    return text[:1].upper() + text[1:].lower()


def _parse_header(line: str, prev_blank: bool = True):
    """Return (kind, title, rest) if ``line`` is a section header, else None.

    Headers sit at indent <= 1 (<= 2 for "Keyword:" after a blank line).
    ALL-CAPS keywords always count; mixed-case ones only with a colon or
    alone on the line ("Example usage:" and "See also a, b" also count).
    Any other line of 2+ ALL-CAPS words (or one ending in ":") is a rubric.
    """
    if not line.strip():
        return None
    indent = len(line) - len(line.lstrip(" "))
    text = line.strip()
    m = _KEYWORD_RE.match(text)
    if m and indent <= 2:
        kw, sep, rest = m.group("kw"), (m.group("sep") or "").strip(), m.group("rest")
        kind, is_caps = _classify_keyword(kw), kw.upper() == kw
        if rest and not sep and not rest[0].isspace():
            m = None  # "Inputsfoo", "Notes," ...
        elif indent == 2 and not (prev_blank and sep == ":"):
            m = None
        elif sep or not rest.strip():
            return kind, _title(kw), rest.strip()
        elif kind == "seealso" and _IDENT_LIST_RE.fullmatch(rest.strip()):
            return kind, _title(kw), rest.strip()
        elif is_caps and re.fullmatch(r"\s+\([^()]*\)\s*", rest):
            return kind, _title(kw), ""  # "OPTIONS (name-value)"
        elif is_caps and indent <= 1 and rest.startswith("  "):
            return kind, _title(kw), rest.strip()
        else:
            m2 = re.fullmatch(r"((?:\s+[A-Za-z]+){1,2})\s*:\s*(.*)", rest)
            if m2:  # "Example usage:"
                return kind, _title(kw + m2.group(1)), m2.group(2).strip()
    if indent > 1:
        return None
    g = _GENERIC_HEADER_RE.match(text)
    if g:
        title, colon, rest = g.group("title"), g.group("colon"), g.group("rest")
        words = [w for w in title.split() if w != "&"]
        if not any(len(w) > 1 for w in words) or (rest.strip() and not colon):
            return None
        if len(words) >= 2 or colon:
            return "rubric", _title(title), rest.strip()
    return None


# ==========================================================================
# Inline text
# ==========================================================================


def escape(text: str) -> str:
    """Escape reST inline markup characters in prose."""
    text = text.replace("\\", "\\\\")
    for ch in "*`|":
        text = text.replace(ch, "\\" + ch)
    text = re.sub(r"(\w)_(?=\W|$)", r"\1\\_", text)  # "word_" is a reference
    text = re.sub(r"(?:(?<=\s)|^)_(?=\w)", r"\\_", text)  # "_word" starts a target
    return text


def _literal(text: str) -> str:
    text = text.strip()
    return escape(text) if (not text or "`" in text) else f"``{text}``"


def _join(lines) -> str:
    return re.sub(r"\s+", " ", " ".join(l.strip() for l in lines)).strip()


# ==========================================================================
# Usage lines
# ==========================================================================


def _is_usage(text: str, names: set[str]) -> bool:
    """True for a call line such as "[a, b] = NAME(x, y)" or "NAME(x)"."""
    if not names:
        return False
    alt = "|".join(re.escape(n) for n in sorted(names, key=len, reverse=True))
    m = re.match(
        rf"^(?P<lhs>(?:\[[^\]]*\]|[\w.]+)\s*=\s*)?(?:\w+\.)*(?:{alt})\b(?P<rest>.*)$",
        text.strip(),
        re.IGNORECASE,
    )
    if not m:
        return False
    rest = m.group("rest").strip()
    if not rest:
        return bool(m.group("lhs"))
    if not rest.startswith("("):
        return False
    depth = 0
    for i, ch in enumerate(rest):
        depth += ch == "("
        depth -= ch == ")"
        if depth == 0:
            after = rest[i + 1 :].strip()
            return after in ("", ";", ",") or after.startswith("%")
    return rest.endswith("...")  # call continued on the next line


def _drop_usage(lines, names):
    """Remove usage lines, including calls continued with "..." on later lines."""
    out, k = [], 0
    while k < len(lines):
        ind, t = lines[k]
        if t.strip():
            # join this line with up to 4 continuation lines and test each prefix
            joined, end = t.rstrip(), k
            hit = _is_usage(joined, names)
            while not hit and joined.endswith("...") and end + 1 < len(lines) and end - k < 4:
                end += 1
                joined = joined[:-3].rstrip() + " " + lines[end][1].strip()
                hit = _is_usage(joined, names)
            if hit:
                # also skip lines the call itself continues onto
                while lines[end][1].rstrip().endswith("...") and end + 1 < len(lines):
                    end += 1
                k = end + 1
                continue
        out.append((ind, t))
        k += 1
    return out


# ==========================================================================
# Parameter / return entries
# ==========================================================================

_NAME = r"(?:'[A-Za-z_]\w*'|[A-Za-z_][\w.]*(?:\([^()]*\))?)"  # 'Name': a name-value option
_NAMES = rf"{_NAME}(?:\s*(?:,|&|\band\b)\s*{_NAME})*"
_ENTRY_RE = re.compile(
    rf"^(?P<names>{_NAMES})(?P<opt>\s*\(optional\))?"
    r"(?:(?P<sep>\s*-\s+|\s*:\s+|\s*:$|\s{2,}|\s+(?=\.[A-Za-z_]))(?P<rest>.*))?$",
    re.IGNORECASE,
)
_EQ_ENTRY_RE = re.compile(  # "Y(:,j) = y(tvec(j))     state at each time"
    r"^(?P<lhs>(?P<name>[A-Za-z_]\w*)(?:\([^()]*\))?)\s*=\s*(?P<rhs>\S.*?)\s{2,}(?P<desc>\S.*)$"
)
_FIELD_RE = re.compile(
    r"^(?P<names>\.[A-Za-z_]\w*(?:\([^)]*\))?(?:\s*,\s*\.[A-Za-z_]\w*(?:\([^)]*\))?)*)"
    r"(?:(?P<sep>\s*-\s+|\s*:\s+|\s+)(?P<rest>.*))?$"
)
_DIMS_RE = re.compile(r"^\s*[\w.+\-*/^ ]+?\s*x\s*[\w.+\-*/^ ]+?(?:\s*x\s*[\w.+\-*/^ ]+?)*\s*$")


@dataclass
class Node:
    """An entry, struct field or text line (``kind``: entry | field | text)."""

    kind: str
    names: list[str]
    column: int = 0
    desc: list[str] = field(default_factory=list)
    children: list["Node"] = field(default_factory=list)
    desc_col: int | None = None
    prefix: str = ""  # reST put before the description (already escaped)
    explicit: bool = False  # entry written "name - desc" or "name: desc"


def _split_names(names: str) -> list[str]:
    return [n.strip().strip("'") for n in re.split(r"\s*(?:,|&|\band\b)\s*", names) if n.strip()]


def _fmt_dims(d: str) -> str:
    return re.sub(r"\s*x\s*", " x ", d.strip())


def _process_desc(lines: list[str]):
    """Join description lines; return (type, escaped text)."""
    text = " ".join(l.strip() for l in lines if l.strip())
    dtype = default = None
    m = re.match(r"^\(([^()]+)\)\s{2,}(.*)$", text)  # "(value)  desc"
    if m:
        default, text = m.group(1).strip(), m.group(2)
    m = re.match(r"^(?:\[([^\]]*)\]|\(([^()]*)\))\s*(.*)$", text)  # leading "[dims]" or "(dims)"
    dims = (m.group(1) if m.group(1) is not None else m.group(2)) if m else None
    if dims is not None and _DIMS_RE.match(dims):
        dtype, text = _fmt_dims(dims), m.group(3)
    last = lines[-1].rstrip() if lines else ""
    m = re.search(r"(?:\s{2,}\[([^\]]+)\]|\s\[([^\s\]]{1,10})\])\s*$", last)  # trailing bracket
    if m:
        val = (m.group(1) or m.group(2)).strip()
        text = re.sub(r"\s*\[" + re.escape(m.group(1) or m.group(2)) + r"\]\s*$", "", text)
        if _DIMS_RE.match(val):
            dtype = dtype or _fmt_dims(val)
        else:
            default = val
    out = escape(re.sub(r"\s+", " ", text).strip())
    if default is not None:
        out = out.rstrip(".").rstrip()
        out = (out + ". " if out else "") + f"Default: {_literal(default)}."
    return dtype, out


def _parse_field_line(text: str, column: int):
    """Parse ".a desc" (or ".sub  .a desc"); return a field Node or None."""
    m = _FIELD_RE.match(text)
    if not m:
        return None
    node = Node("field", [n.strip() for n in m.group("names").split(",")], column)
    rest = m.group("rest") or ""
    if rest.lstrip().startswith(".") and m.group("sep"):
        inner = _parse_field_line(rest.lstrip(), column + len(text) - len(rest.lstrip()))
        if inner:
            node.children.append(inner)
            return node
    if rest.strip():
        node.desc.append(rest)
        node.desc_col = column + len(text) - len(rest.lstrip())
    return node


def _match_entry(stripped: str):
    """Return (names, rest, prefix) if ``stripped`` starts an entry, else None."""
    m = _EQ_ENTRY_RE.match(stripped)
    if m:
        lhs, rhs = m.group("lhs"), m.group("rhs").strip()
        return [m.group("name")], m.group("desc"), f"``{lhs} = {rhs}``: ", True
    m = _ENTRY_RE.match(stripped)
    if not m:
        return None
    alone = stripped == (m.group("names") + (m.group("opt") or "")).strip()
    if m.group("sep") is None and not alone:
        return None
    prefix = "(optional) " if m.group("opt") else ""
    explicit = bool(m.group("sep") and m.group("sep").strip() in ("-", ":"))
    return _split_names(m.group("names")), m.group("rest") or "", prefix, explicit


def _parse_entries(lines):
    """Parse a Parameters/Returns body -> (entries, lead lines, trailing lines).

    Entries start at the indent of the body's first line. Lines before the
    first entry are "lead" prose; less-indented (or same-indent, non-entry)
    lines after it are "trailing" prose shown after the fields.
    """
    first = next(((i, t) for i, t in lines if t.strip()), None)
    if first is None:
        return [], [], []
    entry_indent = first[0]
    entries: list[Node] = []
    stack: list[Node] = []
    lead, trailing = [], []
    prev_blank = True
    for ind, text in lines:
        if not text.strip():
            prev_blank = True
            (trailing if entries else lead).append((ind, ""))
            continue
        was_blank, prev_blank = prev_blank, False
        stripped = text.strip()
        match = _match_entry(stripped)
        # a new entry: at the entry indent or, once a list has started, within
        # 2 columns of it with an explicit "name - desc" / "name: desc"
        at_entry = ind == entry_indent or bool(
            match and match[3] and entries and abs(ind - entry_indent) <= 2
        )
        # "name desc" with one space, aligned with the previous description
        # column, in a list whose entries are separated by spaces only
        prev = entries[-1] if entries else None
        if not match and prev and ind == entry_indent and prev.desc_col and not prev.explicit:
            m1 = re.match(r"^([A-Za-z_]\w*) (\S.*)$", stripped)
            if m1 and ind + len(m1.group(1)) + 1 == prev.desc_col:
                match, at_entry = ([m1.group(1)], m1.group(2), "", False), True
        if match and at_entry:
            names, rest, prefix, explicit = match
            rest_col = ind + len(stripped) - len(rest.lstrip()) if rest.strip() else None
            if len(names) == 1 and "." in names[0] and "(" not in names[0]:
                # "root.sub  .field ..." -> entry "root", field ".sub"
                root, sub = names[0].split(".", 1)
                node = Node("entry", [root], ind)
                fld = Node("field", ["." + sub], ind + len(root))
                node.children.append(fld)
                stack = [fld]
                inner = _parse_field_line(rest.lstrip(), rest_col) if rest.lstrip().startswith(".") else None
                if inner:
                    fld.children.append(inner)
                    stack.append(inner)
                elif rest.strip():
                    fld.desc.append(rest)
                    fld.desc_col = rest_col
                entries.append(node)
                continue
            node = Node("entry", names, ind, prefix=prefix, explicit=explicit)
            entries.append(node)
            stack = []
            if rest.lstrip().startswith(".") and rest_col is not None:
                fld = _parse_field_line(rest.lstrip(), rest_col)
                if fld:
                    node.children.append(fld)
                    stack = [fld] + fld.children[-1:]
                    continue
            if rest.strip():
                node.desc.append(rest)
                node.desc_col = rest_col
            continue
        if not entries:
            lead.append((ind, text))
            continue
        if ind <= entry_indent:  # same-indent prose or outdented: after the fields
            trailing.append((ind, text))
            continue
        entry = entries[-1]
        if stripped.startswith("."):
            fld = _parse_field_line(stripped, ind)
            if fld:
                while stack and stack[-1].column >= ind:
                    stack.pop()
                (stack[-1] if stack else entry).children.append(fld)
                stack.append(fld)
                stack += fld.children[-1:]
                continue
        if stack:
            last = stack[-1]
            if ind > last.column and not was_blank:
                last.desc.append(text)
                continue
            # non-field text among fields -> italic bullet at this level
            while stack and stack[-1].column >= ind:
                stack.pop()
            parent = stack[-1] if stack else entry
            if parent.children and parent.children[-1].kind == "text" and not was_blank:
                parent.children[-1].desc.append(text)
            else:
                parent.children.append(Node("text", [], ind, [text]))
            continue
        if entry.children:
            entry.children.append(Node("text", [], ind, [text]))
            continue
        entry.desc.append(text)
    return entries, lead, trailing


def _render_children(children: list[Node], indent: str) -> list[str]:
    out: list[str] = []
    for ch in children:
        if ch.kind == "text":
            out.append(f"{indent}* *{escape(_join(ch.desc))}*")
        else:
            dtype, desc = _process_desc(ch.desc)
            line = ", ".join(f"``{n}``" for n in ch.names) + (f" (*{escape(dtype)}*)" if dtype else "")
            if desc:
                line += f" -- {desc}"
            out.append(f"{indent}* {line}")
        if ch.children:
            out += [""] + _render_children(ch.children, indent + "  ") + [""]
    return out


def _render_params(entries: list[Node]) -> list[str]:
    out: list[str] = []
    for e in entries:
        name = ",".join(e.names)
        dtype, desc = _process_desc(e.desc)
        out.append(f":param {name}: {e.prefix}{desc}".rstrip())
        if e.children:
            out += [""] + _render_children(e.children, "    ") + [""]
        if dtype:
            out.append(f":type {name}: {escape(dtype)}")
    return out


def _render_returns(entries: list[Node]) -> list[str]:
    if not entries:
        return []
    out = [":returns:", ""]
    for e in entries:
        dtype, desc = _process_desc(e.desc)
        line = f"**{escape(','.join(e.names))}**" + (f" (*{escape(dtype)}*)" if dtype else "")
        if desc or e.prefix:
            line += f" -- {e.prefix}{desc}"
        out.append(f"    * {line}")
        if e.children:
            out += [""] + _render_children(e.children, "      ") + [""]
    return out


# ==========================================================================
# Prose
# ==========================================================================

_MATHY_RE = re.compile(r"=|~|∝|->|>=|<=|\bsum\b")
_CODE_RE = re.compile(r"=\s*[\w.]+\(|;\s*$|^\s*(?:figure|plot|for|if|while|end|disp|load|save)\b")
_BASE_CODE_RE = re.compile(r";\s*(?:%.*)?$|\.\.\.\s*$")  # code at the paragraph's own indent
_BULLET_RE = re.compile(r"^(?:[-*•]|\d+[.)])\s+")


def _code_block(lines, lang: str | None = None) -> list[str]:
    nonblank = [i for i, t in lines if t.strip()]
    if not nonblank:
        return []
    base = min(nonblank)
    content = [(" " * (i - base) + t.strip()) if t.strip() else "" for i, t in lines]
    while content and not content[-1]:
        content.pop()
    while content and not content[0]:
        content.pop(0)
    if lang is None:
        lang = "matlab" if any(_CODE_RE.search(c) for c in content) else "none"
    return [f".. code-block:: {lang}", ""] + ["    " + c if c else "" for c in content] + [""]


def _paragraph(lines: list[str]) -> list[str]:
    txt = _join(lines)
    if not txt:
        return []
    txt = escape(txt)
    if re.match(r"^(?:\.\.|[-+*]\s|\d+[.)]\s|#\.\s|[A-Za-z][.)]\s)", txt):
        txt = "\\" + txt  # not markup in the source
    return [txt, ""]


def _render_list(items: list[list[str]], numbered: bool) -> list[str]:
    return [f"{f'{n}.' if numbered else '*'} {escape(_join(item))}" for n, item in enumerate(items, 1)] + [""]


def render_prose(lines, base: int | None = None) -> list[str]:
    """Render (indent, text) lines as reST paragraphs, lists and code blocks.

    ``base`` is the prose indent; lines indented beyond it are code. It
    defaults to the smallest indent in ``lines``.
    """
    nonblank = [i for i, t in lines if t.strip()]
    if not nonblank:
        return []
    base = min(nonblank) if base is None else min(base, min(nonblank))
    blocks, cur = [], []
    for ind, t in lines:
        if t.strip():
            cur.append((ind, t.strip()))
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)

    out: list[str] = []
    for block in blocks:
        if all(i > base for i, _ in block) and not _BULLET_RE.match(block[0][1]):
            out += _code_block(block)  # fully indented block
            continue
        j, para = 0, []
        while j < len(block):
            ind, t = block[j]
            if _BULLET_RE.match(t):
                out += _paragraph(para)
                para, items, col = [], [], ind
                numbered = t[0].isdigit()
                while j < len(block):
                    ind2, t2 = block[j]
                    m2 = _BULLET_RE.match(t2)
                    if m2 and ind2 == col:
                        items.append([t2[m2.end():]])
                    elif ind2 > col and items:
                        items[-1].append(t2)
                    else:
                        break
                    j += 1
                out += _render_list(items, numbered)
                continue
            if ind <= base and _BASE_CODE_RE.search(t):
                # code written at the paragraph's own indent ("x = f(y);")
                out += _paragraph(para)
                para, run = [], []
                while j < len(block) and (
                    _BASE_CODE_RE.search(block[j][1]) or block[j][0] > base
                    or (run and _BASE_CODE_RE.search(run[-1][1]) and run[-1][1].rstrip().endswith("..."))
                ):
                    run.append(block[j])
                    j += 1
                out += _code_block(run, "matlab")
                continue
            if ind > base and j > 0:
                run, k = [], j
                while k < len(block) and block[k][0] > base and not _BULLET_RE.match(block[k][1]):
                    run.append(block[k])
                    k += 1
                prev = block[j - 1][1].rstrip()
                listing = any(re.match(rf"^{_NAME}\s+-\s+\S", r[1]) for r in run)  # "x - desc" lines
                if prev.endswith(":") or listing or _MATHY_RE.search(prev) or any(_MATHY_RE.search(r[1]) for r in run):
                    out += _paragraph(para)
                    para = []
                    out += _code_block(run)
                else:
                    para += [r[1] for r in run]  # hanging indent: same paragraph
                j = k
                continue
            para.append(t)
            j += 1
        out += _paragraph(para)
    return out


# ==========================================================================
# Whole help block
# ==========================================================================


def _seealso(lines, known: set[str]) -> list[str]:
    refs = []
    for w in re.split(r"[\s,;]+", " ".join(t for _, t in lines)):
        w = w.strip().rstrip(".")
        if not w:
            continue
        parts = w.split(".")
        short = parts[-1]
        lower = {k.lower(): k for k in known}
        if short in known:
            refs.append(f":mat:obj:`{short}`")
        elif short.lower() in lower:
            refs.append(f":mat:obj:`{lower[short.lower()]}`")
        elif len(parts) >= 2 and parts[-2] in known:  # Class.method
            refs.append(f":mat:obj:`{parts[-2]}.{short} <{parts[-2]}.{short}>`")
        elif re.fullmatch(r"[\w.]+", w):
            refs.append(f"``{w}``")
        else:
            refs.append(escape(w))
    return [".. rubric:: See also", "", ", ".join(refs), ""] if refs else []


def _summary(tagged, short: str):
    """Split off the summary -> (summary lines, rest, name aliases, kept subject)."""
    aliases = {short.lower()} if short else set()
    first = tagged[0][1].strip()
    if _parse_header(tagged[0][1]):
        return [], tagged, aliases, ""
    m = re.match(r"^(?P<fn>function\s+)?(?P<tok>(?:\w+\.)*\w+)(?P<sep>\s*:\s*|\s{2,}|\s+|$)(?P<rest>.*)$", first)
    if not m:
        return [], tagged, aliases, ""
    tok, sep, rest = m.group("tok").split(".")[-1], m.group("sep"), m.group("rest").strip()
    caps = tok.isupper() and len(tok) > 2
    loud = ":" in sep or len(sep) >= 2 or not sep
    is_name = (short and tok.lower() == short.lower()) or (caps and (loud or m.group("fn")))
    if not is_name:
        return [], tagged, aliases, ""
    aliases.add(tok.lower())
    if _is_usage(rest, aliases) or re.match(r"^\[.*\]\s*=", rest):
        rest = ""
    subject = ""
    if rest and re.match(r"^(?:is|are|was|takes|provides|returns)\b", rest) and not loud:
        subject = short or tok  # "NAME is a class ..." keeps its subject
    summary = [rest] if rest else []
    k = 1
    while k < len(tagged) and tagged[k][1].strip():
        t = tagged[k][1]
        if _parse_header(t, False) or _is_usage(t, aliases):
            break
        summary.append(t.strip())
        k += 1
    return summary, tagged[k:], aliases, subject


def _author_text(lines: list[str]) -> str:
    out = ""
    for l in lines:
        if not out:
            out = l
        elif out.rstrip().endswith((":", ",", "&", "and")):
            out += " " + l
        else:
            out += "; " + l
    return out


def convert(raw: list[str], name: str = "", known: set[str] | None = None) -> list[str]:
    """Convert MATLAB help text (list of lines, ``%`` stripped) to reST lines."""
    known = known or set()
    short = name.split(".")[-1]
    lines = [l.expandtabs(4).rstrip() for l in raw]
    lines = [l for l in lines if not _BANNER_RE.match(l)]
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        return []
    base = min(len(l) - len(l.lstrip()) for l in lines if l.strip())
    tagged = [(len(l) - len(l.lstrip()) - base, l[base:]) if l.strip() else (0, "") for l in lines]

    summary, rest_lines, aliases, subject = _summary(tagged, short)

    # ---- sections
    sections = [["prose", "", []]]
    prev_blank = True
    for ind, t in rest_lines:
        h = _parse_header(t, prev_blank) if t.strip() else None
        prev_blank = not t.strip()
        if h:
            kind, title, rest = h
            sections.append([kind, title, [(None, rest)] if rest else []])
            continue
        sections[-1][2].append((ind, t))
    for sec in sections:
        body = sec[2]
        nb = [i for i, t in body if t.strip() and i is not None]
        first_ind = min(nb) if nb else 0
        body = [(first_ind if i is None else i, t) for i, t in body]
        if sec[0] != "example":
            body = _drop_usage(body, aliases)
        sec[2] = body

    body_indents = [i for _k, _t, body in sections for i, t in body if t.strip()]
    body_base = min(body_indents) if body_indents else 0

    # ---- render
    out = _paragraph(summary)
    if out and subject:
        out[0] = f"``{subject}`` {out[0]}"
    params, returns = [], []
    returns_note = ""
    pre_fields, after_fields, tail, author = [], [], [], []
    started_tail = False
    for kind, title, body in sections:
        has_text = any(t.strip() for _, t in body)
        if kind in ("prose", "usage"):
            # the description before the first header has no header to set its
            # indent: use the smallest indent of the whole body, so that a
            # fully indented block (equation, diagram) stays code
            (tail if started_tail else out).extend(render_prose(body, None if started_tail else body_base))
        elif kind in ("params", "returns", "properties"):
            text = " ".join(t.strip() for _, t in body if t.strip())
            if not has_text or re.match(r"^none\b", text, re.IGNORECASE):
                if kind == "returns" and has_text:  # "none (draws into the axes)"
                    returns_note = "None" + escape(text[4:]).rstrip(".") + "."
                started_tail = True
                continue
            if kind == "properties":
                entries, lead, trailing = _parse_entries(body)
                tail += [f".. rubric:: {escape(title)}", ""] + render_prose(lead)
                tail += [l[4:] if l.startswith("    ") else l for l in _render_returns(entries)[2:]] + [""]
                tail += render_prose(trailing)
                started_tail = True
                continue
            entries, lead, trailing = _parse_entries(body)
            if not entries:
                tail += [f".. rubric:: {escape(title)}", ""] + render_prose(lead + trailing)
            else:
                (params if kind == "params" else returns).extend(entries)
                pre_fields += render_prose(lead)
                after_fields += render_prose(trailing)
            started_tail = True
        elif kind == "example":
            if has_text:
                tail += [f".. rubric:: {escape(title)}", ""] + _code_block(body, "matlab")
            started_tail = True
        elif kind == "seealso":
            tail += _seealso(body, known)
            started_tail = True
        elif kind == "author":
            author += [t.strip() for _, t in body if t.strip()]
            started_tail = True
        else:  # NOTES, REFERENCES, any other header
            tail += [f".. rubric:: {escape(title)}", ""] + render_prose(body)
            started_tail = True

    out += pre_fields
    fields = _render_params(params) + _render_returns(returns)
    if returns_note and not returns:
        fields.append(f":returns: {returns_note}")
    if fields:
        if out and out[-1]:
            out.append("")
        out += fields + [""]
    out += after_fields + tail
    if author:
        out += [f"*Author:* {escape(_author_text(author))}", ""]
    while out and not out[-1]:
        out.pop()
    return out


# ==========================================================================
# Sphinx glue
# ==========================================================================

_DUMP: list[dict] = []


def _known_names(app) -> set[str]:
    names = set()
    for root, _dirs, files in os.walk(app.config.matlab_src_dir or ".", followlinks=True):
        names.update(f[:-2] for f in files if f.endswith(".m"))
    return names


def _verbatim(lines: list[str]) -> list[str]:
    body = [l.rstrip() for l in lines]
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    return [".. code-block:: none", ""] + ["    " + l for l in body] + [""] if body else []


def process_docstring(app, what, name, obj, options, lines):
    raw = list(lines)
    if not any(l.strip() for l in raw):
        return
    known = getattr(app, "_matlab_help_known", None)
    if known is None:
        known = app._matlab_help_known = _known_names(app)
    try:
        new, fallback = convert(raw, name, known), False
    except Exception as exc:  # noqa: BLE001 - any parse failure -> verbatim
        logger.warning("matlab_help: could not parse the help text of %s (%s); shown verbatim", name, exc)
        new, fallback = _verbatim(raw), True
    lines[:] = new + [""]
    if os.environ.get("MATLAB_HELP_DUMP"):
        _DUMP.append({"what": what, "name": name, "raw": raw, "rst": new, "fallback": fallback})


def _write_dump(app, exc):
    path = os.environ.get("MATLAB_HELP_DUMP")
    if path:
        Path(path).write_text(json.dumps(_DUMP, indent=1))


def setup(app):
    app.connect("autodoc-process-docstring", process_docstring)
    app.connect("build-finished", _write_dump)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
