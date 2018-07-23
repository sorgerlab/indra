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
        return

    def make_report(self):
        doc = SimpleDocTemplate(self.name + '.pdf', pagesize=letter,
                                rightMargin=72, leftMargin=72,
                                topMargin=72, bottomMargin=18)
        doc.build(self.story)
        return doc

    def add_story_text(self, text, style='Normal', space=None, fontsize=12):
        if space is None:
            space=(1,12)
        ptext = '<fond size=%d>%s</font>' % (fontsize, text)
        self.story.append(Paragraph(ptext, self.styles[style]))
        self.story.append(Spacer(*space))
        return

    def add_story_image(self, image_path, width=None, height=None):
        if width is not None:
            width = width*inch
        if height is not None:
            height = height*inch
        im = Image(image_path, width, height)
        self.story.append(im)
