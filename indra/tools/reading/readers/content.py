import zlib
import shutil

from os import path


class Content(object):
    """An object to regularize the content passed to the readers.

    To use this class, use one of the two constructor methods:
     - `from_file` : use content from a file on the filesystem.
     - `from_string` : Pass a string (or bytes) directly as content.

    This class also regularizes the handling of id's and formats, as well as
    allowing for decompression and decoding, in the manner standard in the INDRA
    project.
    """
    def __init__(self, id, format, compressed=False, encoded=False):
        self.file_exists = False
        self.compressed = compressed
        self.encoded = encoded
        self._id = id
        self._format = format
        self._text = None
        self._fname = None
        self._location = None
        self._raw_content = None
        return

    def __repr__(self):
        return ('Content(id=\'%s\', path=\'%s\')'
                % (self.get_id(), self.get_filepath()))

    @classmethod
    def from_file(cls, file_path, compressed=False, encoded=False):
        """Create a content object from a file path."""
        file_id = '.'.join(path.basename(file_path).split('.')[:-1])
        file_format = file_path.split('.')[-1]
        content = cls(file_id, file_format, compressed, encoded)
        content.file_exists = True
        content._location = path.dirname(file_path)
        return content

    @classmethod
    def from_string(cls, id, format, raw_content, compressed=False,
                    encoded=False):
        """Create a Content object from string/bytes content."""
        content = cls(id, format, compressed, encoded)
        content._raw_content = raw_content
        return content

    def _load_raw_content(self):
        if self.file_exists and self._raw_content is None:
            with open(self.get_filepath(), 'r') as f:
                self._raw_content = f.read()
        return

    def change_id(self, new_id):
        """Change the id of this content."""
        self._load_raw_content()
        self._id = new_id
        self.get_filename(renew=True)
        self.get_filepath(renew=True)
        return

    def change_format(self, new_format):
        """Change the format label of this content.

        Note that this does NOT actually alter the format of the content, only
        the label.
        """
        self._load_raw_content()
        self._format = new_format
        self.get_filename(renew=True)
        self.get_filepath(renew=True)
        return

    def set_location(self, new_location):
        """Set/change the location of this content.

        Note that this does NOT change the actual location of the file. To do
        so, use the `copy_to` method.
        """
        self._load_raw_content()
        self._location = new_location
        self.get_filepath(renew=True)
        return

    def is_format(self, *formats):
        """Check the format of this content."""
        return any([self._format == fmt for fmt in formats])

    def get_id(self):
        return self._id

    def get_format(self):
        return self._format

    def get_text(self):
        """Get the loaded, decompressed, and decoded text of this content."""
        self._load_raw_content()
        if self._text is None:
            assert self._raw_content is not None
            ret_cont = self._raw_content
            if self.compressed:
                ret_cont = zlib.decompress(ret_cont, zlib.MAX_WBITS+16)
            if self.encoded:
                ret_cont = ret_cont.decode('utf-8')
            self._text = ret_cont
        assert self._text is not None
        return self._text

    def get_filename(self, renew=False):
        """Get the filename of this content.

        If the file name doesn't already exist, we created it as {id}.{format}.
        """
        if self._fname is None or renew:
            self._fname = '%s.%s' % (self._id, self._format)
        return self._fname

    def get_filepath(self, renew=False):
        """Get the file path, joining the name and location for this file.

        If no location is given, it is assumed to be "here", e.g. ".".
        """
        if self._location is None or renew:
            self._location = '.'
        return path.join(self._location, self.get_filename())

    def copy_to(self, location, fname=None):
        if fname is None:
            fname = self.get_filename()
        fpath = path.join(location, fname)
        if self.file_exists and not self._raw_content:
            shutil.copy(self.get_filepath(), fpath)
        else:
            with open(fpath, 'w') as f:
                f.write(self.get_text())
        self._fname = fname
        self._location = location
        self.file_exists = True
        return fpath
