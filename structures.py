import copy


class Vec:

  def __init__(self, contents=None):
    """constructor defaults to empty vector
           accepts list of elements to initialize a vector object with the
           given list
        """
    if contents is None:
      contents = []
    self.elements = contents
    return

  def __abs__(self):
    """Overloads the built-in function abs(v)
            returns the Euclidean norm of vector v
        """
    return sum([e**2 for e in self.elements])**0.5

  def __add__(self, other):
    """Overloads the + operation to support Vec + Vec
         raises ValueError if vectors are not same length
        """
    if type(other) != Vec:
      raise ValueError(f"Vec + {type(other)} is not defined.")
    if len(self.elements) == len(other.elements):
      n = len(self.elements)
      return Vec([self.elements[i] + other.elements[i] for i in range(n)])
    else:
      raise ValueError("ERROR: Vectors must be same length")

  def __mul__(self, other):
    """Overloads the * operator to support
            - Vec * Vec (dot product) raises ValueError if vectors are not same length in the case of dot product
            - Vec * float (component-wise product)
            - Vec * int (component-wise product)

        """
    if type(other) == Vec:  # define dot product
      if len(self.elements) == len(other.elements):
        n = len(self.elements)
        return sum([self.elements[i] * other.elements[i] for i in range(n)])
      else:
        raise ValueError("ERROR: Vectors must be same length")
    elif type(other) == float or type(other) == int or type(
        other) == complex:  # scalar-vector multiplication
      return Vec([x * other for x in self.elements])
    else:
      raise ValueError(f"Vec * {type(other)} is not supported.")

  def __rmul__(self, other):
    """Overloads the * operation to support
            - float * Vec
            - int * Vec
        """
    if type(other) == float or type(other) == int or type(other) == complex:
      return Vec([x * other for x in self.elements])
    else:
      raise ValueError(
          f"ERROR: {type(other)} * {type(self)} is not supported.")

  def __truediv__(self, other):
    if type(other) == complex or type(other) == int or type(other) == float:
      return Vec([x / other for x in self.elements])
    else:
      raise ValueError(f"Vec / {type(other)} is not defined.")

  def __str__(self):
    """returns string representation of this Vec object"""
    return str(self.elements)  # does NOT need further implementation

  def __sub__(self, other):
    if len(self.elements) == len(other.elements):
      n = len(self.elements)
      return Vec([self.elements[i] - other.elements[i] for i in range(n)])
    else:
      raise ValueError("ERROR: Vectors must be same length")

  def __getitem__(self, i):
    return self.elements[i]

  def __eq__(self, other):
    return self.elements == other.elements

  def norm(self, p):
    return sum([(abs(x)**p) for x in self.elements])**(1 / p)

  def dim(self):
    return len(self.elements)


class Matrix:

  def __init__(self, rows=[]):
    self.rowsp = rows
    self.colsp = []
    self._set_colsp()
    return

  def dim(self):
    m = len(self.rowsp)
    n = len(self.colsp)
    return (m, n)

  def _set_colsp(self):
    """HELPER METHOD: Resets the column space according to the existing row space"""
    self.colsp = []
    n = len(self.rowsp[0])
    m = len(self.rowsp)
    for j in range(n):
      col = []
      for i in range(m):
        col.append(self.rowsp[i][j])
      self.colsp.append(col)
    return

  def _set_rowsp(self):
    """HELPER METHOD: Resets the row space according to the existing column space"""
    self.rowsp = []
    n = len(self.colsp)
    m = len(self.colsp[0])
    for i in range(m):
      row = []
      for j in range(n):
        row.append(self.colsp[j][i])
      self.rowsp.append(row)

  def transpose(self):
    return Matrix(copy.deepcopy(self.colsp))

  def __str__(self):
    """returns string representation of this Matrix object"""
    return str(self.rowsp)  # does NOT need further implementation

  #FIXME: COPY AND PASTE THE REST OF YOUR MATRIX METHODS HERE
