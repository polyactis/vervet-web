from vervetserver.tests import *

class TestAccessionController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='Accession', action='index'))
        # Test response...
