import logging

from pylons import request, response, session, tmpl_context as c, url
from pylons.controllers.util import abort, redirect

from vervetserver.lib.base import BaseController, render

log = logging.getLogger(__name__)

class AccessionController(BaseController):

    def index(self):
        # Return a rendered template
        #return render('/Accession.mako')
        # or, return a string
        return 'Hello World'