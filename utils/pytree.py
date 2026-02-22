#!/usr/bin/env python3
import os
import argparse
import fnmatch

def parse_args():
    parser = argparse.ArgumentParser(
        description="Minimal tree-like directory viewer (no install required)"
    )
    parser.add_argument(
        "path", nargs="?", default=".", help="Root path (default: .)"
    )
    parser.add_argument(
        "-L", "--level", type=int, default=None, help="Max display depth"
    )
    parser.add_argument(
        "-d", "--dirs-only", action="store_true", help="List directories only"
    )
    parser.add_argument(
        "-a", "--all", action="store_true", help="Show hidden files"
    )
    parser.add_argument(
        "-I", "--ignore", action="append", default=[], help="Ignore pattern"
    )
    return parser.parse_args()


def ignored(name, patterns):
    return any(fnmatch.fnmatch(name, p) for p in patterns)


def walk(path, prefix="", depth=0, args=None):
    if args.level is not None and depth >= args.level:
        return

    try:
        entries = os.listdir(path)
    except PermissionError:
        return

    entries = sorted(entries)

    filtered = []
    for e in entries:
        if not args.all and e.startswith("."):
            continue
        if ignored(e, args.ignore):
            continue
        full = os.path.join(path, e)
        if args.dirs_only and not os.path.isdir(full):
            continue
        filtered.append(e)

    for i, name in enumerate(filtered):
        full = os.path.join(path, name)
        last = (i == len(filtered) - 1)
        connector = "└── " if last else "├── "
        print(prefix + connector + name)

        if os.path.isdir(full):
            extension = "    " if last else "│   "
            walk(full, prefix + extension, depth + 1, args)


def main():
    args = parse_args()
    root = os.path.abspath(args.path)

    print(os.path.basename(root) or root)
    walk(root, args=args)


if __name__ == "__main__":
    main()
