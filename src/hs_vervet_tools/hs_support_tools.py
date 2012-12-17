# module with some support tools

class hsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def safeWrite(self,fname,write_str,exist_action='a'):
    """
    exist_action...
    'w'... overwrite existing file
    'a'... append to existing file
    'p'... pass if the file exists
    """
    import os
    if exist_action=='p':
        if os.path.isfile(fname):
            pass
        else:
            with open(fname, 'w') as f:
                f.write(write_str)
    elif exist_action=='a':
        with open(fname, 'a') as f:
            f.write(write_str)
    elif exist_action=='w':
        with open(fname, 'w') as f:
            f.write(write_str)
    else:
        print "variable existAction must be w (overwrite), a (append), or p (pass)"