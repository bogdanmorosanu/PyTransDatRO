import math

DEG_TO_RAD_FACTOR = math.pi / 180

def sexa_to_rad(s):
    """Returns the value in radians of a sexagesimal representation of an angle
    The format of the representation is: [DD MM SS.SSS]

    :param s: DMS angle
    :type s: string

    .. warnings also:: This function doesn't make any check on the input. 
        It is designed to be used internaly, so valid in put is expected
    """
    sgn = -1 if s[0] == '-' else 1
    values = s.split(' ')
    d = abs(float(values[0]))
    m = float(values[1])
    s = float(values[2])
    return  sgn * DEG_TO_RAD_FACTOR * (s / 3600 + m / 60 + d)

def rad_to_sexa(r):
    """[summary]

    :param r: angle in radians
    :type r: float
    """
    sgn = '-' if r < 0 else ''
    tmp_d = abs(r / DEG_TO_RAD_FACTOR)
    d = int(tmp_d)
    tmp_reminder = tmp_d - d
    m = int(tmp_reminder * 60)
    tmp_reminder = (tmp_reminder * 60) - m
    s = tmp_reminder * 60
    return f'{sgn}{str(d)} {str(m)} {str(s)}'
    



def ppm_to_unity(ppm):
    """Converts a ppm value to its unity value equivalent
    
    :param ppm: ppm scale value
    :type ppm: float   

    :return: the unity scale value of the input ppm value
    :rtype: float
    """
    return 1 + ppm * 1E-6



