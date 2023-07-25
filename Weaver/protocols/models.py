from django.db import models
import uuid
from inventory.custom.general import COLORS
from organization.models import Project
from django.contrib.auth.models import User

REACTIVE_STATE = (
    (0, 'Solid'),
    (1, 'Liquid')
)

CONCENTRATION_UNITS = (
    ('mol', 'M'),
    ('mmol', 'mM'),
    ('grlt', 'mg / ml (= g/l)'),
    ('vvp', 'Volume / Volume %'),
    ('wvp', 'Weight / Volume %'),
)


class Reactive(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(User, on_delete=models.CASCADE, editable=False)
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
    owner = models.ForeignKey(User, on_delete=models.CASCADE, editable=False)
    reactive = models.ForeignKey(Reactive, on_delete=models.CASCADE)
    concentration = models.FloatField()
    concentration_unit = models.CharField(choices=CONCENTRATION_UNITS, max_length=4)

    class Meta:
        ordering = ["reactive"]

    def __str__(self):
        return self.reactive.name + " - " + self.concentration.__str__() + " " + self.concentration_units

    @property
    def concentration_units(self):
        return dict(CONCENTRATION_UNITS).get(self.concentration_unit.__str__())


class TableFilter(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=20, blank=True)
    color = models.CharField(choices=COLORS, max_length=10, blank=True, null=True)

    def __str__(self):
        return self.name


class Recipe(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(User, on_delete=models.CASCADE, editable=False)
    name = models.CharField(max_length=200)
    components = models.ManyToManyField(Component, blank=True)
    ph = models.CharField(max_length=10, blank=True, null=True)
    description = models.CharField(max_length=1000, blank=True, null=True)
    category = models.ManyToManyField(TableFilter, blank=True, help_text="Use CTRL for multiple select")
    shared_to_project = models.ManyToManyField(Project, blank=True, help_text="Use CTRL for multiple select")

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['name']
