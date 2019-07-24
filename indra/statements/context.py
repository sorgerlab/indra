import datetime


__all__ = ['Context', 'BioContext', 'WorldContext', 'RefContext', 'TimeContext',
           'MovementContext']


class Context(object):
    """An abstract class for Contexts."""
    @classmethod
    def from_json(cls, jd):
        context_type = jd.get('type')
        if context_type == 'bio':
            return BioContext.from_json(jd)
        elif context_type == 'world':
            return WorldContext.from_json(jd)
        elif context_type == 'movement':
            return MovementContext.from_json(jd)
        else:
            raise ValueError('Unknown context type %s' % context_type)


class BioContext(Context):
    """An object representing the context of a Statement in biology.

    Parameters
    ----------
    location : Optional[RefContext]
        Cellular location, typically a sub-cellular compartment.
    cell_line : Optional[RefContext]
        Cell line context, e.g., a specific cell line, like BT20.
    cell_type : Optional[RefContext]
        Cell type context, broader than a cell line, like macrophage.
    organ : Optional[RefContext]
        Organ context.
    disease : Optional[RefContext]
        Disease context.
    species : Optional[RefContext]
        Species context.
    """
    def __init__(self, location=None, cell_line=None, cell_type=None,
                 organ=None, disease=None, species=None):
        self.location = location
        self.cell_line = cell_line
        self.cell_type = cell_type
        self.organ = organ
        self.disease = disease
        self.species = species

    attrs = ['location', 'cell_line', 'cell_type', 'organ', 'disease',
             'species']

    def __eq__(self, other):
        return all([getattr(self, attr, None) == getattr(other, attr, None)
                    for attr in self.attrs])

    def __ne__(self, other):
        return not self.__eq__(other)

    def __bool__(self):
        return any([getattr(self, attr, None) is not None
                    for attr in self.attrs])

    def __nonzero__(self):
        return self.__bool__()

    @classmethod
    def from_json(cls, jd):
        # For all the attributes, we deserialize them if they have a value,
        # and make a dict that can be passed to the constructor
        ref_contexts = {attr: (RefContext.from_json(jd.get(attr))
                               if jd.get(attr) else None)
                        for attr in cls.attrs}
        bs = cls(**ref_contexts)
        return bs

    def to_json(self):
        jd = {attr: getattr(self, attr).to_json() for attr in self.attrs
              if getattr(self, attr, None) is not None}
        jd['type'] = 'bio'
        return jd

    def __str__(self):
        pieces = ['%s=%s' % (attr, getattr(self, attr)) for attr in self.attrs
                  if getattr(self, attr, None) is not None]
        args = ', '.join(pieces)
        return 'BioContext(%s)' % args

    def __repr__(self):
        return self.__str__()


class WorldContext(Context):
    """An object representing the context of a Statement in time and space.

    Parameters
    ----------
    time : Optional[TimeContext]
        A TimeContext object representing the temporal context of the
        Statement.
    geo_location : Optional[RefContext]
        The geographical location context represented as a RefContext
    """
    def __init__(self, time=None, geo_location=None):
        self.time = time
        self.geo_location = geo_location

    def __eq__(self, other):
        return (self.time == other.time) and \
               (self.geo_location == other.geo_location)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __bool__(self):
        return self.time is not None or self.geo_location is not None

    def __nonzero__(self):
        return self.__bool__()

    @classmethod
    def from_json(cls, jd):
        time_entry = jd.get('time')
        time = TimeContext.from_json(time_entry) if time_entry else None
        geo_entry = jd.get('geo_location')
        geo_location = RefContext.from_json(geo_entry) if geo_entry else None
        return cls(time=time, geo_location=geo_location)

    def to_json(self):
        jd = {'type': 'world',
              'time': self.time.to_json() if self.time else None,
              'geo_location': (self.geo_location.to_json()
                               if self.geo_location else None)}
        return jd

    def __str__(self):
        pieces = []
        if self.time:
            pieces.append('time=%s' % self.time)
        if self.geo_location:
            pieces.append('geo_location=%s' % self.geo_location)
        args = ', '.join(pieces)
        return 'WorldContext(%s)' % args

    def __repr__(self):
        return self.__str__()


class RefContext(object):
    """An object representing a context with a name and references.

    Parameters
    ----------
    name : Optional[str]
        The name of the given context. In some cases a text name will not be
        available so this is an optional parameter with the default being
        None.
    db_refs : Optional[dict]
        A dictionary where each key is a namespace and each value is an
        identifier in that namespace, similar to the db_refs associated with
        Concepts/Agents.
    """
    def __init__(self, name=None, db_refs=None):
        self.name = name
        self.db_refs = {} if db_refs is None else db_refs

    def __eq__(self, other):
        return (self.name == other.name) and \
               (self.db_refs == other.db_refs)

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_json(self):
        jd = {'name': self.name, 'db_refs': self.db_refs}
        return jd

    @classmethod
    def from_json(cls, jd):
        rc = cls(name=jd.get('name'), db_refs=jd.get('db_refs'))
        return rc

    def __str__(self):
        pieces = []
        if self.name:
            pieces.append('name="%s"' % self.name)
        if self.db_refs:
            pieces.append('db_refs=%s' % self.db_refs)
        args = ', '.join(pieces)
        return 'RefContext(%s)' % args

    def __repr__(self):
        return self.__str__()


class TimeContext(object):
    """An object representing the time context of a Statement

    Parameters
    ----------
    text : Optional[str]
        A string representation of the time constraint, typically as seen in
        text.
    start : Optional[datetime]
        A `datetime` object representing the start time
    end : Optional[datetime]
        A `datetime` object representing the end time
    duration : int
        The duration of the time constraint in seconds
    """
    def __init__(self, text=None, start=None, end=None, duration=None):
        self.text = text
        self.start = start
        self.end = end
        self.duration = duration

    def refinement_of(self, other):
        # If both starts are given
        if self.start and other.start:
            # If this started earlier, it can't be a refinement
            if self.start <= other.start:
                return False
            # If it ended later, it can't be a refinement
            elif self.end and other.end:
                if self.end >= other.end:
                    return False
            # If the other end is not given, it's a refinement
            elif self.end and not other.end:
                return True
            # If neither end is given, it's also a refinement at this point
            return True
        # If we have a start time and other doesn't, it's a refinement
        elif self.start and not other.start:
            return True
        # Otherwise it's not a refinement
        return False

    def __eq__(self, other):
        return (self.text == other.text) and \
               (self.start == other.start) and \
               (self.end == other.end) and \
               (self.duration == other.duration)

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_json(self):
        def date_to_str(date):
            if date is None:
                return None
            else:
                return date.strftime('%Y-%m-%dT%H:%M')
        jd = {'text': self.text,
              'start': date_to_str(self.start),
              'end': date_to_str(self.end)}
        if self.duration is not None:
            jd['duration'] = self.duration
        return jd

    @classmethod
    def from_json(cls, jd):
        def date_from_str(date_str):
            if not date_str:
                return None
            try:
                dt = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M')
            except Exception as e:
                return None
            return dt
        tc = cls(text=jd.get('text'),
                 start=date_from_str(jd.get('start')),
                 end=date_from_str(jd.get('end')),
                 duration=jd.get('duration'))
        return tc

    def __str__(self):
        pieces = []
        if self.text:
            pieces.append('text="%s"' % self.text)
        if self.start:
            pieces.append('start=%s' % self.start)
        if self.end:
            pieces.append('end=%s' % self.end)
        if self.duration:
            pieces.append('duration=%s' % self.duration)
        args = ', '.join(pieces)
        return 'TimeContext(%s)' % args

    def __repr__(self):
        return self.__str__()


class MovementContext(Context):
    """An object representing the context of a movement between start and end
    points in time.

    Parameters
    ----------
    locations : Optional[list[dict]
        A list of dictionaries each containing a RefContext object representing
        geographical location context and its role (e.g. 'origin',
        'destination', etc.)
    time : Optional[TimeContext]
        A TimeContext object representing the temporal context of the
        Statement.
    """
    def __init__(self, locations=None, time=None):
        self.locations = locations if locations else []
        self.time = time

    def __eq__(self, other):
        return (self.locations == other.locations) and \
               (self.time == other.time)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __bool__(self):
        return self.locations is not None or self.time is not None

    def __nonzero__(self):
        return self.__bool__()

    @classmethod
    def from_json(cls, jd):
        locations_entry = jd.get('locations')
        if locations_entry:
            locations = [
                {'location': RefContext.from_json(location['location']),
                 'role': location['role']} for location in locations_entry]
        else:
            locations = None
        time_entry = jd.get('time')
        if time_entry:
            time = TimeContext.from_json(time_entry)
        else:
            time = None
        return cls(locations=locations, time=time)

    def to_json(self):
        if self.locations:
            locations_json = [
                {'location': location['location'].to_json(),
                 'role': location['role']} for location in self.locations]
        else:
            locations_json = []
        jd = {'type': 'movement',
              'locations': locations_json,
              'time': self.time.to_json() if self.time else None}
        return jd

    def __str__(self):
        pieces = []
        if self.locations:
            locs = []
            for location in self.locations:
                locs.append('%s(%s)' % (location['location'], location['role']))
            locations = ', '.join(locs)
            pieces.append('locations=%s' % locations)
        if self.time:
            pieces.append('time=%s' % self.time)
        args = ', '.join(pieces)
        return 'MovementContext(%s)' % args

    def __repr__(self):
        return self.__str__()
