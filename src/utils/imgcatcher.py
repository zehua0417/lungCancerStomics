import os
from matplotlib.backends.backend_pdf import PdfPages


class ImgCatcher:
    """
    Class to save images in a pdf file
    """
    def __init__(self, filename):
        self.filename = filename
        self.pdf = PdfPages(filename)

    def save(self):
        # Check if the directory exists, if not, create it
        try:
            directory = os.path.dirname(self.filename)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)
            self.pdf.savefig()
        except Exception as e:
            print(f"Error: {e}")

    def close(self):
        self.pdf.close()
        del self

    def save_and_close(self):
        self.save()
        self.close()
