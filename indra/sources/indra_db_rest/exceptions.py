class IndraDBRestClientError(Exception):
    pass


class IndraDBRestResponseError(IndraDBRestClientError):
    pass


class IndraDBRestAPIError(IndraDBRestClientError):
    def __init__(self, resp):
        self.status_code = resp.status_code
        if hasattr(resp, 'text'):
            self.reason = resp.text
        else:
            self.reason = resp.reason
        reason_lines = self.reason.splitlines()
        if len(reason_lines) > 1:
            ws = '\n'
            for i, line in enumerate(reason_lines):
                reason_lines[i] = '  ' + line
        else:
            ws = ' '
        fmtd_reason = '\n'.join(reason_lines)
        self.resp = resp
        Exception.__init__(self, ('Got bad return code %d:%s%s'
                                  % (self.status_code, ws, fmtd_reason)))
        return
