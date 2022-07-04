# Weaver | Open plasmid admin
Open source web administrator for plasmids + strains with tools for molecular biology.

Django + React dev

## Demo

Visit our demo at: [weaver-demo.meristem.cl/](https://weaver-demo.meristem.cl/)

Default super user credentials:
|User|Password  |
|--|--|
| root | pass |


## Installing
### Pre-requisites
Make sure you are running python > 3.6. To check your python version run:

`python -v`

Install following packages via python package installer (pip):

`pip3 install django-mathfilters django-shortuuidfield django-multiselectfield Bio plotly pandas`

### Clone repo
Clone repo:

`git clone https://github.com/meristem-group/Weaver`

## Running
### Start web server

Enter project folder:

`cd Weaver/Weaver`

Start the server:

`python3 manage.py runserver`

Access local server : [http://127.0.0.1:8000/](http://127.0.0.1:8000/)

Default super user credentials:
|User|Password  |
|--|--|
| root | pass |

### Creating objects
Groups & Users from AUTHENTICATION AND AUTHORIZATION and Boxs, Locations & Strains from INVENTORY can only be edited from [admin site](http://127.0.0.1:8000/admin/). All other elements can be created and edited from front-end.