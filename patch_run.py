#!/usr/bin/env python3
"""Guarded legacy patch helper for source/run.cc."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUN_SOURCE = ROOT / "source" / "run.cc"


OLD_LOOP_PATTERN = (
    r"        LA::MPI::Vector damage_prev_iter\(locally_owned_dofs_damage(.*?)"
    r"        completely_distributed_solution_damage_old =\n"
    r"\s*completely_distributed_solution_damage;"
)

NEW_LOOP = """        bool stoppingCriterion = false;
        unsigned int iteration = 0;

        while (stoppingCriterion == false)
        {
            assemble_system_elastic();
            solve_linear_system_elastic();
            locally_relevant_solution_elastic.update_ghost_values();

            assemble_system_damage();
            solve_linear_system_damage();

            locally_relevant_solution_damage.update_ghost_values();

            if (iteration > 0)
                stoppingCriterion = check_convergence ();

            completely_distributed_solution_elastic_old =
                locally_relevant_solution_elastic;
            completely_distributed_solution_damage_old =
                locally_relevant_solution_damage;

            if (stoppingCriterion == false)
            {
                refine_grid(boundary_value);
            }
            iteration = iteration + 1;
        } // while loop for converged solution ends

        update_history_field();"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply the legacy run.cc patch. Omit this flag for a dry check.",
    )
    args = parser.parse_args()

    content = RUN_SOURCE.read_text(encoding="utf-8")
    patched, count = re.subn(OLD_LOOP_PATTERN, NEW_LOOP, content, flags=re.DOTALL)

    if count == 0:
        print(f"No legacy loop pattern found in {RUN_SOURCE}")
        return 0

    if not args.apply:
        print(f"Dry run: legacy loop pattern found {count} time(s); no files changed.")
        return 0

    RUN_SOURCE.write_text(patched, encoding="utf-8")
    print(f"Patched {RUN_SOURCE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
