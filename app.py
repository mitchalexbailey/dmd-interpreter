import os
from urllib.parse import quote

from flask import Flask, render_template

from interpreter import views

app = Flask(
    __name__,
    template_folder=os.path.join(os.path.dirname(__file__), "interpreter", "templates"),
    static_folder=os.path.join(os.path.dirname(__file__), "interpreter", "static"),
)

app.add_url_rule("/", view_func=views.index, methods=["GET"])
app.add_url_rule("/results", view_func=views.results, methods=["POST"])
app.add_url_rule("/api/interpret", view_func=views.api_interpret, methods=["GET", "POST"])

app.jinja_env.filters["urlencode"] = lambda value: quote(value or "")


@app.errorhandler(404)
def not_found(_error):
    return render_template("404.html"), 404


application = app

if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8000"))
    app.run(host="0.0.0.0", port=port, debug=True)
