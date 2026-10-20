"""The command line guide says what is required, what defaults, and what is off.

A guide that has drifted from the parser is worse than none, so the facts in it are
read back from the parser here.
"""

from __future__ import annotations

import re
from pathlib import Path

from promoter_tf_drugs import cli

GUIDE = (Path(__file__).resolve().parents[1] / "docs" / "cli-guide.md"
         ).read_text(encoding="utf-8")

# Options the guide leaves to the command reference on purpose.
LEFT_TO_THE_REFERENCE = {
    "--matrix-id", "--promoter-fasta", "--lincs-file", "--signature-file",
    "--cache-dir", "--max-requests", "--activity-permutations", "--help",
    "--condition", "--tissue",
}


def _options():
    return [a for a in cli.build_parser()._actions
            if a.option_strings and a.option_strings[0] != "-h"]


def test_the_only_required_option_is_the_one_the_guide_names():
    required = [a.option_strings[0] for a in _options()
                if getattr(a, "required", False)]
    assert required == ["--gene"]
    assert "Only `--gene` is required" in GUIDE


def test_the_defaults_table_matches_the_parser():
    block = GUIDE.split("<!-- defaults:start -->")[1].split(
        "<!-- defaults:end -->")[0]
    rows = re.findall(r"^\| `(--[\w-]+)` \| `([^`]*)` \|", block, re.M)
    assert len(rows) >= 14
    parser_defaults = {a.option_strings[0]: a.default for a in _options()}
    for opt, shown in rows:
        assert opt in parser_defaults, opt
        assert str(parser_defaults[opt]) == shown, (
            f"{opt}: the guide says {shown}, the parser says "
            f"{parser_defaults[opt]}")


def test_every_switch_that_turns_something_off_is_in_the_guide():
    missing = []
    for dest in cli._NEGATED_ARGS:
        flag = "--" + dest.replace("_", "-")
        if flag not in GUIDE:
            missing.append(flag)
    assert not missing, missing


def test_every_switch_that_is_off_by_default_is_in_the_guide():
    neg = set(cli._NEGATED_ARGS)
    missing = [a.option_strings[0] for a in _options()
               if a.nargs == 0 and a.dest not in neg
               and a.option_strings[0] not in LEFT_TO_THE_REFERENCE
               and a.option_strings[0] not in GUIDE]
    assert not missing, f"add {missing} to docs/cli-guide.md"


def test_every_option_that_takes_a_file_or_a_name_is_in_the_guide():
    missing = [a.option_strings[0] for a in _options()
               if a.nargs != 0 and a.default is None
               and a.option_strings[0] not in LEFT_TO_THE_REFERENCE
               and a.option_strings[0] not in GUIDE]
    assert not missing, f"add {missing} to docs/cli-guide.md"


def test_the_refusals_the_guide_lists_are_real(capsys):
    import pytest
    for argv in (["--gene", "X", "--tf-activity"],
                 ["--gene", "X", "--activity-permutations", "5"],
                 ["--gene", "X", "--scenario", "toxic_gain", "--direction",
                  "up"]):
        with pytest.raises(SystemExit):
            cli.main(argv)
    capsys.readouterr()
