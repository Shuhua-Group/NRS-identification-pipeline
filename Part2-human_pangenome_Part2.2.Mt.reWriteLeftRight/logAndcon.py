
class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module
    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        print(msg, file=self.log_fh)
        print(msg)

class Config(object):
    def __init__(self, fh):
        self.con_fh = open(fh, 'r')
        self.configs = self._readConfig()

    def _readConfig(self):
        configs = {}
        for line in self.con_fh.readlines():
            line = line.strip()
            if line[0] != '#':
                line = line.split('=')
                #line = [''.join(s.split()) for s in line ]
                line = [s.strip() for s in line]
                configs[line[0]] = line[1]
        return configs
    def getValue(self, value):
        return self.configs.get(value)