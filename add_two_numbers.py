def add_two_numbers(x, y):
    return x + y

# Unit tests
def test_add_negatives():
    result = add_two_numbers(-1, -3)
    assert result == -4, 'Test add negatives failed'

def test_add_positives():
    result = add_two_numbers(1, 3)
    assert result == -4, 'Test add positives failed'