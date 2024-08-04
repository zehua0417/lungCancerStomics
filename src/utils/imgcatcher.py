from matplotlib.backends.backend_pdf import PdfPages

class ImgCatcher:
    def __init__(self, filename):
        self.filename = filename
        self.pdf = PdfPages(filename)

    def save(self):
        try:
            self.pdf.savefig()
        except Exception as e:
            print(f"Error: {e}")

    def close(self):
        self.pdf.close()
        # del self
    
    def save_and_close(self):
        self.save()
        self.close()