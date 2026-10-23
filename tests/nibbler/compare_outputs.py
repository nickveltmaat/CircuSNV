#!/usr/bin/env python3
"""
Compare the output directories of two CircuSNV runs at record level.

Usage:
    compare_outputs.py <output_dir_a> <output_dir_b> [options]

Every file that exists under either directory is matched by its relative path and
compared as a multiset of records, so record order does not matter:
  - VCF (.vcf, .vcf.gz): header and data lines;
  - text files: lines;
  - Excel (.xlsx): rows, per sheet (needs openpyxl, which the pipeline venv has);
  - other binary files: SHA-256 of the content.

Volatile lines are ignored and counted: VCF '##' header lines whose key mentions a
date or command (fileDate, *Command, GATKCommandLine, ...) and header-like lines or
Excel rows that contain a date, a time, an absolute path or the word 'command'.
Index files, OpenCRAVAT databases and log files are skipped (see --exclude).

The report contains counts only. Sample directory names are masked as S1, S2, ...
unless --no-mask is given, so the report can be shared without identifiers.

Exit status: 0 = no differences, 1 = differences found, 2 = usage or read error.
"""

import argparse
import collections
import fnmatch
import gzip
import hashlib
import os
import re
import sys

DEFAULT_EXCLUDES = ['*.csi', '*.tbi', '*.bai', '*.crai', '*.sqlite', '*.log', '*.err', '*.png']

CALLER_RE = re.compile(r'_(VarDict|Mutect2|LoFreq|SiNVICT)(?:_RAW)?\.vcf')
VOLATILE_KEY_RE = re.compile(r'date|command|cmdline', re.IGNORECASE)
DATE_RE = re.compile(r'\b(?:19|20)\d{2}-\d{2}-\d{2}\b|\b\d{2}:\d{2}:\d{2}\b')
PATH_RE = re.compile(r'(?:^|[\s=:"\'(,</])/[^\s/"\',>]+/')
COMMAND_RE = re.compile(r'\bcommand\b', re.IGNORECASE)

TABLE_COLUMNS = ['files', 'same', 'differ', 'only-A', 'only-B', 'skipped',
                 'records-A', 'records-B', 'rec-only-A', 'rec-only-B', 'ignored']


class ReadError(Exception):
    """A file could not be read. The message never contains the path, which may be internal."""


def reason(exc):
    return getattr(exc, 'strerror', None) or type(exc).__name__


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Compare two CircuSNV output directories at record level.')
    parser.add_argument('dir_a', help='output directory of version A (e.g. the baseline)')
    parser.add_argument('dir_b', help='output directory of version B (e.g. the candidate)')
    parser.add_argument('--exclude', action='append', default=[], metavar='GLOB',
                        help='skip files whose relative path or name matches GLOB (repeatable)')
    parser.add_argument('--no-mask', action='store_true',
                        help='show real sample directory names (do not share that report)')
    parser.add_argument('--mapping-file', metavar='PATH',
                        help='write the sample-to-mask mapping to PATH (keep it on the cluster)')
    parser.add_argument('--show-records', type=int, default=0, metavar='N',
                        help='print up to N differing records per file; these contain variants, '
                             'so never share that part of the report')
    parser.add_argument('--skip-xlsx', action='store_true',
                        help='skip Excel files instead of failing when openpyxl is missing')
    parser.add_argument('--show-paths', action='store_true',
                        help='print the two directory paths in the report header')
    return parser.parse_args(argv)


def list_files(root):
    files = set()
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in filenames:
            full = os.path.join(dirpath, name)
            files.add(os.path.relpath(full, root).replace(os.sep, '/'))
    return files


def top_level_dirs(root):
    return {name for name in os.listdir(root) if os.path.isdir(os.path.join(root, name))}


class Masker(object):
    """Replaces sample directory names by S1, S2, ... in paths and records."""

    def __init__(self, samples, enabled):
        self.enabled = enabled
        self.mapping = {name: 'S%d' % (i + 1) for i, name in enumerate(sorted(samples))}
        # Longest names first, so a name that is a prefix of another is not replaced inside it.
        self.order = sorted(self.mapping, key=len, reverse=True)

    def __call__(self, text):
        if not self.enabled:
            return text
        for name in self.order:
            text = text.replace(name, self.mapping[name])
        return text


def category_of(rel):
    parts = rel.split('/')
    if len(parts) == 1:
        return 'top-level'
    if len(parts) == 2:
        match = CALLER_RE.search(parts[1])
        return 'caller:' + match.group(1) if match else 'sample-root'
    return parts[1]


def is_excluded(rel, patterns):
    name = rel.rsplit('/', 1)[-1]
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(name, p) for p in patterns)


def file_kind(path):
    lower = path.lower()
    if lower.endswith('.xlsx'):
        return 'xlsx'
    if lower.endswith('.vcf') or lower.endswith('.vcf.gz'):
        return 'vcf'
    try:
        if lower.endswith('.gz'):
            with gzip.open(path, 'rb') as fh:
                chunk = fh.read(8192)
        else:
            with open(path, 'rb') as fh:
                chunk = fh.read(8192)
    except (OSError, EOFError) as exc:
        raise ReadError(reason(exc))
    return 'binary' if b'\0' in chunk else 'text'


def is_volatile_line(line):
    if line.startswith('##'):
        key = line[2:].split('=', 1)[0]
        if VOLATILE_KEY_RE.search(key):
            return True
    return bool(DATE_RE.search(line) or PATH_RE.search(line) or COMMAND_RE.search(line))


def open_text(path):
    if path.lower().endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace', newline='')
    return open(path, 'r', encoding='utf-8', errors='replace', newline='')


def read_text_records(path, kind):
    """Return (Counter of lines, number of data lines, number of ignored lines)."""
    records = collections.Counter()
    data = 0
    ignored = 0
    try:
        with open_text(path) as fh:
            for raw in fh:
                line = raw.rstrip('\r\n')
                if line.startswith('#'):
                    if is_volatile_line(line):
                        ignored += 1
                        continue
                    records[line] += 1
                else:
                    records[line] += 1
                    data += 1
    except (OSError, EOFError) as exc:
        raise ReadError(reason(exc))
    return records, data, ignored


def cell_text(value):
    return '' if value is None else str(value)


def read_xlsx_records(path):
    """Return (Counter of 'sheet<TAB>cells' rows, number of data rows, number of ignored rows)."""
    import openpyxl  # imported here so VCF and text comparison works without it

    records = collections.Counter()
    data = 0
    ignored = 0
    try:
        workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except Exception as exc:  # openpyxl raises several unrelated exception types
        raise ReadError(reason(exc))
    try:
        for sheet in workbook.worksheets:
            for row in sheet.iter_rows(values_only=True):
                cells = [cell_text(v) for v in row]
                while cells and cells[-1] == '':
                    cells.pop()
                if not cells:
                    continue
                joined = '\t'.join(cells)
                if is_volatile_line(joined):
                    ignored += 1
                    continue
                records[sheet.title + '\t' + joined] += 1
                data += 1
    finally:
        workbook.close()
    return records, data, ignored


def sha256_of(path):
    digest = hashlib.sha256()
    try:
        with open(path, 'rb') as fh:
            for block in iter(lambda: fh.read(1 << 20), b''):
                digest.update(block)
    except OSError as exc:
        raise ReadError(reason(exc))
    return digest.hexdigest()


def compare_file(path_a, path_b, have_openpyxl, skip_xlsx):
    """Return a dict describing the comparison, or None when the file type is skipped."""
    kind = file_kind(path_a)
    if kind == 'xlsx':
        if not have_openpyxl:
            if skip_xlsx:
                return None
            raise ReadError('openpyxl is not available (run with the pipeline venv, '
                            'or pass --skip-xlsx)')
        rec_a, data_a, ign_a = read_xlsx_records(path_a)
        rec_b, data_b, ign_b = read_xlsx_records(path_b)
    elif kind == 'binary':
        same = sha256_of(path_a) == sha256_of(path_b)
        return {'same': same, 'records_a': 0, 'records_b': 0, 'only_a': collections.Counter(),
                'only_b': collections.Counter(), 'ignored': 0, 'binary': True}
    else:
        rec_a, data_a, ign_a = read_text_records(path_a, kind)
        rec_b, data_b, ign_b = read_text_records(path_b, kind)
    only_a = rec_a - rec_b
    only_b = rec_b - rec_a
    return {'same': not only_a and not only_b, 'records_a': data_a, 'records_b': data_b,
            'only_a': only_a, 'only_b': only_b, 'ignored': ign_a + ign_b, 'binary': False}


def empty_row():
    return collections.OrderedDict((col, 0) for col in TABLE_COLUMNS)


def print_table(stats):
    width = max([len('category')] + [len(c) for c in stats])
    header = 'category'.ljust(width) + ''.join(col.rjust(12) for col in TABLE_COLUMNS)
    print(header)
    print('-' * len(header))
    totals = empty_row()
    for category in sorted(stats):
        row = stats[category]
        print(category.ljust(width) + ''.join(str(row[col]).rjust(12) for col in TABLE_COLUMNS))
        for col in TABLE_COLUMNS:
            totals[col] += row[col]
    print('-' * len(header))
    print('total'.ljust(width) + ''.join(str(totals[col]).rjust(12) for col in TABLE_COLUMNS))


def main(argv):
    args = parse_args(argv)
    for directory in (args.dir_a, args.dir_b):
        if not os.path.isdir(directory):
            sys.stderr.write('ERROR: not a directory: %s\n' % directory)
            return 2

    try:
        import openpyxl  # noqa: F401
        have_openpyxl = True
    except ImportError:
        have_openpyxl = False

    excludes = DEFAULT_EXCLUDES + args.exclude
    mask = Masker(top_level_dirs(args.dir_a) | top_level_dirs(args.dir_b), not args.no_mask)
    if args.mapping_file:
        with open(args.mapping_file, 'w') as fh:
            for name in sorted(mask.mapping):
                fh.write('%s\t%s\n' % (mask.mapping[name], name))

    files_a = list_files(args.dir_a)
    files_b = list_files(args.dir_b)
    stats = collections.defaultdict(empty_row)
    differing = []
    only_in_a = []
    only_in_b = []
    errors = []

    for rel in sorted(files_a | files_b):
        category = category_of(rel)
        row = stats[category]
        row['files'] += 1
        if is_excluded(rel, excludes):
            row['skipped'] += 1
            continue
        if rel not in files_b:
            row['only-A'] += 1
            only_in_a.append(rel)
            continue
        if rel not in files_a:
            row['only-B'] += 1
            only_in_b.append(rel)
            continue
        try:
            result = compare_file(os.path.join(args.dir_a, rel), os.path.join(args.dir_b, rel),
                                  have_openpyxl, args.skip_xlsx)
        except ReadError as exc:
            errors.append('%s: %s' % (mask(rel), exc))
            continue
        if result is None:
            row['skipped'] += 1
            continue
        n_only_a = sum(result['only_a'].values())
        n_only_b = sum(result['only_b'].values())
        row['records-A'] += result['records_a']
        row['records-B'] += result['records_b']
        row['rec-only-A'] += n_only_a
        row['rec-only-B'] += n_only_b
        row['ignored'] += result['ignored']
        if result['same']:
            row['same'] += 1
        else:
            row['differ'] += 1
            differing.append((rel, result, n_only_a, n_only_b))

    print('CircuSNV output comparison')
    if args.show_paths:
        print('A: %s' % args.dir_a)
        print('B: %s' % args.dir_b)
    print('Sample directories: %d%s' % (len(mask.mapping), ' (masked)' if mask.enabled else ''))
    print('Skipped patterns: %s' % ' '.join(excludes))
    print('')
    print_table(stats)

    if differing:
        print('\nDiffering files:')
        for rel, result, n_only_a, n_only_b in differing:
            if result['binary']:
                print('  %s: binary content differs' % mask(rel))
            else:
                print('  %s: %d record(s) only in A, %d only in B'
                      % (mask(rel), n_only_a, n_only_b))
            if args.show_records > 0 and not result['binary']:
                for label, counter in (('A', result['only_a']), ('B', result['only_b'])):
                    for record in sorted(counter)[:args.show_records]:
                        print('      only in %s: %s' % (label, mask(record)))
    if only_in_a:
        print('\nFiles only in A:')
        for rel in only_in_a:
            print('  %s' % mask(rel))
    if only_in_b:
        print('\nFiles only in B:')
        for rel in only_in_b:
            print('  %s' % mask(rel))
    if errors:
        print('\nErrors:')
        for message in errors:
            print('  %s' % message)

    if errors:
        print('\nResult: ERROR (%d file(s) could not be compared)' % len(errors))
        return 2
    if differing or only_in_a or only_in_b:
        print('\nResult: DIFFERENT')
        return 1
    print('\nResult: IDENTICAL')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
