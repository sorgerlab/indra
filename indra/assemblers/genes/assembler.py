# -*- coding: utf-8 -*-

"""An assembler for gene lists and related gene-centric outputs."""

from pathlib import Path
from typing import List, Optional, Set, Union

from indra.statements import Statement

__all__ = [
    "GeneAssembler",
]


class GeneAssembler:
    """The Gene assembler assembles INDRA Statements into a gene list.

    This graph can then be used with the GSEA software from the Broad
    Institute, or output as a visualization with Ideogram.

    Parameters
    ----------
    stmts :
        A list of INDRA Statements to be added to the assembler's list
        of Statements.

    Example
    -------
    To create an ideogram HTML file, use the following code:

    .. code-block:: python

        stmts = ...
        assembler = GeneAssembler(stmts)
        assembler.make_model()
        assembler.save_model('ideogram.html', format='ideogram')
    """

    stmts: List[Statement]
    #: A set of strings representing gene names in the model.
    genes: Set[str]

    def __init__(self, stmts: Optional[List[Statement]] = None):
        self.stmts = stmts or []
        self.genes = set()

    def make_model(self):
        """Assemble the graph from the assembler's list of INDRA Statements."""
        for statement in self.stmts:
            try:
                agents = statement.agent_list()
            except Exception:
                continue
            else:
                for agent in agents:
                    if agent is None:
                        continue
                    if "HGNC" in agent.db_refs:
                        self.genes.add(agent.name)

        return self.genes

    def save_model(self, path: Union[str, Path], format="gsea"):
        """Save the assembled model's string into a file.

        Parameters
        ----------
        path :
            The name of the file to save the SIF into.
        """
        with open(path, "w") as file:
            if format == "gsea":
                print(self.to_gsea_str(), file=file)
            elif format == "ideogram":
                print(self.to_ideogram_html_str(), file=file)
            else:
                raise ValueError(f"Unsupported format: {format}")

    def to_gsea_str(self, first: Optional[str] = "# INDRA Genes") -> str:
        """Get a gene list string of the assembled model."""
        return first + "\n" + "\n".join(sorted(self.genes))

    def to_ideogram_html_str(self) -> str:
        """Get an Ideogram HTML document string of the assembled model."""
        import pydeogram

        return pydeogram.to_html_str(self.genes)
