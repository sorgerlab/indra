from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch


class Reporter(object):
    """A class to produce a simple pdf report using the reportlab package.

    Parameters
    ----------
    name : str
        A name that will be the name of this file, and the title by default.
    """
    def __init__(self, name):
        self.styles = getSampleStyleSheet()
        self.styles.add(ParagraphStyle(name="Justify", alignment=TA_JUSTIFY))
        self.story = []
        self.name = name
        self.title = name
        self.sections = {}
        self.section_headings = []
        return

    def set_title(self, title):
        """Overwrite the default title with a string title of your choosing."""
        self.title = title
        return

    def add_section(self, section_name):
        """Create a section of the report, to be headed by section_name

        Text and images can be added by using the `section` argument of the
        `add_text` and `add_image` methods. Sections can also be ordered by
        using the `set_section_order` method.

        By default, text and images that have no section will be placed after
        all the sections, in the order they were added. This behavior may be
        altered using the `sections_first` attribute of the `make_report`
        method.
        """
        self.section_headings.append(section_name)
        if section_name in self.sections:
            raise ValueError("Section %s already exists." % section_name)
        self.sections[section_name] = []
        return

    def set_section_order(self, section_name_list):
        """Set the order of the sections, which are by default unorderd.

        Any unlisted sections that exist will be placed at the end of the
        document in no particular order.
        """
        self.section_headings = section_name_list[:]
        for section_name in self.sections.keys():
            if section_name not in section_name_list:
                self.section_headings.append(section_name)
        return

    def add_text(self, text, *args, **kwargs):
        """Add text to the document.

        Text is shown on the final document in the order it is added, either
        within the given section or as part of the un-sectioned list of content.

        Parameters
        ----------
        text : str
            The text to be added.
        style : str
            Choose the style of the text. Options include 'Normal', 'Code',
            'Title', 'h1'. For others, see `getSampleStyleSheet` from
            `reportlab.lib.styles`.
        space : tuple (num spaces, font size)
            The number and size of spaces to follow this section of text.
            Default is (1, 12).
        fontsize : int
            The integer font size of the text (e.g. 12 for 12 point font).
            Default is 12.
        alignment : str
            The alignment of the text. Options include 'left', 'right', and
            'center'. Default is 'left'.
        section : str
            (This must be a keyword) Select a section in which to place this
            text. Default is None, in which case the text will be simply be
            added to a default list of text and images.
        """
        # Pull down some kwargs.
        section_name = kwargs.pop('section', None)

        # Actually do the formatting.
        para, sp = self._preformat_text(text, *args, **kwargs)

        # Select the appropriate list to update
        if section_name is None:
            relevant_list = self.story
        else:
            relevant_list = self.sections[section_name]

        # Add the new content to list.
        relevant_list.append(para)
        relevant_list.append(sp)
        return

    def add_image(self, image_path, width=None, height=None, section=None):
        """Add an image to the document.

        Images are shown on the final document in the order they are added,
        either within the given section or as part of the un-sectioned list of
        content.

        Parameters
        ----------
        image_path : str
            A path to the image on the local file system.
        width : int or float
            The width of the image in the document in inches.
        height : int or float
            The height of the image in the document in incehs.
        section : str
            (This must be a keyword) Select a section in which to place this
            image. Default is None, in which case the image will be simply be
            added to a default list of text and images.
        """
        if width is not None:
            width = width*inch
        if height is not None:
            height = height*inch
        im = Image(image_path, width, height)
        if section is None:
            self.story.append(im)
        else:
            self.sections[section].append(im)
        return

    def make_report(self, sections_first=True, section_header_params=None):
        """Create the pdf document with name `self.name + '.pdf'`.

        Parameters
        ----------
        sections_first : bool
            If True (default), text and images with sections are presented first
            and un-sectioned content is appended afterword. If False, sectioned
            text and images will be placed before the sections.
        section_header_params : dict or None
            Optionally overwrite/extend the default formatting for the section
            headers. Default is None.
        """
        full_story = list(self._preformat_text(self.title, style='Title',
                                               fontsize=18, alignment='center'))

        # Set the default section header parameters
        if section_header_params is None:
            section_header_params = {'style': 'h1', 'fontsize': 14,
                                     'alignment': 'center'}

        # Merge the sections and the rest of the story.
        if sections_first:
            full_story += self._make_sections(**section_header_params)
            full_story += self.story
        else:
            full_story += self.story
            full_story += self._make_sections(**section_header_params)

        fname = self.name + '.pdf'
        doc = SimpleDocTemplate(fname, pagesize=letter,
                                rightMargin=72, leftMargin=72,
                                topMargin=72, bottomMargin=18)
        doc.build(full_story)
        return fname

    def _make_sections(self, **section_hdr_params):
        """Flatten the sections into a single story list."""
        sect_story = []
        if not self.section_headings and len(self.sections):
            self.section_headings = self.sections.keys()

        for section_name in self.section_headings:
            section_story = self.sections[section_name]
            line = '-'*20
            section_head_text = '%s %s %s' % (line, section_name, line)
            title, title_sp = self._preformat_text(section_head_text,
                                                   **section_hdr_params)
            sect_story += [title, title_sp] + section_story
        return sect_story

    def _preformat_text(self, text, style='Normal', space=None, fontsize=12,
                        alignment='left'):
        """Format the text for addition to a story list."""
        if space is None:
            space=(1,12)
        ptext = ('<para alignment=\"%s\"><font size=%d>%s</font></para>'
                 % (alignment, fontsize, text))
        para = Paragraph(ptext, self.styles[style])
        sp = Spacer(*space)
        return para, sp
