from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch


class Reporter(object):
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
        self.title = title
        return

    def set_section_order(self, section_name_list):
        self.section_headings = section_name_list[:]
        for section_name in self.sections.keys():
            if section_name not in section_name_list:
                section_name_list.append(section_name)
        return

    def add_section(self, section_name):
        self.section_headings.append(section_name)
        if section_name in self.sections:
            raise ValueError("Section %s already exists." % section_name)
        self.sections[section_name] = []

    def make_report(self, sections_first=True, section_header_params=None):
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
        sect_story = []
        for section_name, section_story in self.sections.items():
            line = '-'*20
            section_head_text = '%s %s %s' % (line, section_name, line)
            title, title_sp = self._preformat_text(section_head_text,
                                                   **section_hdr_params)
            sect_story += [title, title_sp] + section_story
        return sect_story

    def _preformat_text(self, text, style='Normal', space=None, fontsize=12,
                        alignment='left'):
        if space is None:
            space=(1,12)
        ptext = ('<para alignment=\"%s\"><font size=%d>%s</font></para>'
                 % (alignment, fontsize, text))
        para = Paragraph(ptext, self.styles[style])
        sp = Spacer(*space)
        return para, sp

    def add_story_text(self, text, *args, **kwargs):
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

    def add_story_image(self, image_path, width=None, height=None,
                        section=None):
        if width is not None:
            width = width*inch
        if height is not None:
            height = height*inch
        im = Image(image_path, width, height)
        if section is None:
            self.story.append(im)
        else:
            self.sections[section].append(im)
