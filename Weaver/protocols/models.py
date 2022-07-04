from django.db import models
import uuid
from multiselectfield import MultiSelectField

# Create your models here.

REACTIVE_STATE = (
    (0, 'Solid'),
    (1, 'Liquid')
)

CONCENTRATION_UNITS = (
    ('mol', 'M'),
    ('mmol', 'mM'),
    ('grlt', 'gr / lt'),
    ('vvp', 'Volume / Volume %'),
    ('wvp', 'Weight / Volume %'),
)

RECIPE_CATEGORIES = (
    (0, 'General'),
    (1, 'DNA'),
    (2, 'RNA'),
    (3, 'Protein'),
    (4, 'Cell'),
    (5, 'Base editors'),
    (6, 'Protoplasts'),
)


class Reactive(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200)
    state = models.IntegerField(choices=REACTIVE_STATE)
    mm = models.FloatField(blank=True, null=True, help_text="g/mol")
    concentration = models.FloatField(blank=True, null=True)
    concentration_unit = models.CharField(choices=CONCENTRATION_UNITS, max_length=4, blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']


class Component(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    reactive = models.ForeignKey(Reactive, on_delete=models.CASCADE)
    concentration = models.FloatField()
    concentration_unit = models.CharField(choices=CONCENTRATION_UNITS, max_length=4)

    class Meta:
        ordering = ["reactive"]

    def __str__(self):
        return self.reactive.name + " @ " + self.concentration.__str__() + " [" + dict(CONCENTRATION_UNITS).get(
            self.concentration_unit.__str__()) + "]"


class Recipe(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    name = models.CharField(max_length=200)
    components = models.ManyToManyField(Component, blank=True)
    ph = models.CharField(max_length=10, blank=True, null=True)
    description = models.CharField(max_length=1000, blank=True, null=True)
    category = models.IntegerField(choices=RECIPE_CATEGORIES, blank=True, null=True)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']
