# -*- coding: utf-8 -*-

"""Command line interface for DGI-DB."""

from collections import Counter

from tabulate import tabulate

from .processor import DGIProcessor


def main():
    processor = DGIProcessor()
    statements = processor.extract_statements()

    print(f"Number skipped: {processor.skipped}\n")
    print(
        tabulate(
            Counter(
                statement.__class__.__name__ for statement in statements
            ).most_common(),
            headers=["Statement Type", "Count"],
        )
    )

    print()
    print(
        tabulate(
            Counter(
                evidence.annotations["source"]
                for statement in statements
                for evidence in statement.evidence
            ).most_common(),
            headers=["Source", "Count"],
        )
    )

    print()
    print(
        tabulate(
            Counter(
                interaction
                for statement in statements
                for evidence in statement.evidence
                for interaction in evidence.annotations["interactions"].split(",")
            ).most_common(),
            headers=["Interaction", "Count"],
        )
    )


if __name__ == "__main__":
    main()
