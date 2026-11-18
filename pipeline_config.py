#!/usr/bin/env python3
"""
pipeline_config.py
==================
Single source of the site-specific paths used by every script in this
repository (data roots, study spreadsheet, output folders, external tools).

The paths live in config.yaml next to this file (copy config.example.yaml and
edit it). The file is flat YAML: one `key: value` per line, `#` starts a
comment, and a value may refer to another key as ${key}. Relative paths are
resolved against the folder that holds the config file. Another file can be
selected with the environment variable GLYSRC7T_CONFIG.

Only paths are configured here. Analysis parameters (thresholds, seeds,
margins, draw and permutation counts, FDR settings) stay in the scripts and
their command-line flags, exactly as they were run.

Python:   from pipeline_config import cfg;  root = cfg("data_root")
Shell:    python3 pipeline_config.py data_root      (prints one value)
          python3 pipeline_config.py                (prints all values)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ENV_VAR = "GLYSRC7T_CONFIG"
_REF = re.compile(r"\$\{([A-Za-z0-9_]+)\}")
_cache = {}


def config_path():
    """The config file in use: $GLYSRC7T_CONFIG, else config.yaml, else the example."""
    env = os.environ.get(ENV_VAR)
    if env:
        if not os.path.isfile(env):
            sys.exit(f"[config] {ENV_VAR}={env} does not exist")
        return env
    own = os.path.join(HERE, "config.yaml")
    if os.path.isfile(own):
        return own
    example = os.path.join(HERE, "config.example.yaml")
    print(f"[config] config.yaml not found; using the placeholders in {example}. "
          "Copy it to config.yaml and set your paths.", file=sys.stderr)
    return example


def _parse(path):
    raw = {}
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            line = line.split(" #", 1)[0].strip() if not line.lstrip().startswith("#") else ""
            if not line:
                continue
            if ":" not in line:
                raise ValueError(f"{path}:{n}: expected 'key: value', got {line!r}")
            key, val = (s.strip() for s in line.split(":", 1))
            raw[key] = val.strip("'\"")
    return raw


def _resolve(raw, key, base, seen=()):
    if key not in raw:
        raise KeyError(f"'{key}' is not set in {base}")
    if key in seen:
        raise ValueError(f"circular ${{}} reference: {' -> '.join(seen + (key,))}")
    val = _REF.sub(lambda m: _resolve(raw, m.group(1), base, seen + (key,)), raw[key])
    val = os.path.expanduser(val)
    if val and not os.path.isabs(val):
        val = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(base)), val))
    return val


def load():
    """All settings as {key: resolved path}."""
    path = config_path()
    if path not in _cache:
        raw = _parse(path)
        _cache[path] = {k: _resolve(raw, k, path) for k in raw}
    return _cache[path]


def cfg(key):
    """One resolved setting; a missing key is an error that names the config file."""
    settings = load()
    if key not in settings:
        sys.exit(f"[config] '{key}' is not set in {config_path()}")
    return settings[key]


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in ("-h", "--help"):
        print(__doc__.strip())
    elif len(sys.argv) == 2:
        print(cfg(sys.argv[1]))
    elif len(sys.argv) == 1:
        for k, v in load().items():
            print(f"{k}: {v}")
    else:
        sys.exit("usage: python3 pipeline_config.py [key]")
