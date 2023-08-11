def Settings( **kwargs ):
  return {
    'flags': ['-x', 'c++', '-Wall', '-pedantic',
    '-std=c++17', '-I/usr/include', '-I/usr/include/irrlicht',
    '-I/usr/include/bullet'],
  }

