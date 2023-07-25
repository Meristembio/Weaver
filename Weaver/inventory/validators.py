from django.core.exceptions import ValidationError
def clustal_validate(value):
    if value.size > 1048576:
        raise ValidationError("The maximum file size that can be uploaded is 1 MB")
    if not value.name.endswith('.clustal'):
        raise ValidationError("File type must be Clustal (.clustal)")
    else:
        return value
