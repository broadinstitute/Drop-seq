#!/usr/bin/env python3
# MIT License
#
# Copyright 2025 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Sends emails
"""
import copy
import os
import smtplib
import socket
from abc import ABC, abstractmethod
from email.message import EmailMessage
from typing import Optional

try:
    from .models import HtmlEmailMessage, SmtpSettings
except ImportError:
    from models import HtmlEmailMessage, SmtpSettings


class EmailClient(ABC):
    """
    An email client.
    """

    def __init__(self, default_from: Optional[str] = None):
        self._default_from_address: Optional[str] = default_from

    def default_from_address(self) -> str:
        """
        Get the default from email address.
        """
        if not self._default_from_address:
            self._default_from_address = f"{self._default_from_user()}@{self._default_from_domain()}"
        return self._default_from_address

    @abstractmethod
    def send_email(self, message: HtmlEmailMessage) -> None:
        """
        Send an email.
        """
        pass

    def _default_from_user(self):
        """
        Get the default user for the email address.
        """
        return os.getenv("USER")

    def _default_from_domain(self):
        """
        Get the default domain for the email address.
        """
        hostname = socket.gethostname()
        domain = socket.getfqdn()
        if domain.startswith(hostname + "."):
            domain = domain[len(hostname) + 1:]
        return domain


class SmtpEmailClient(EmailClient):
    """
    Send an email via an SMTP server.
    """

    def __init__(self, smtp_settings: SmtpSettings, default_from: Optional[str] = None):
        super().__init__(default_from)
        self.smtp_settings = copy.copy(smtp_settings)

    def send_email(self, message: HtmlEmailMessage) -> None:
        email_message = EmailMessage()
        email_message["From"] = message.email_from
        email_message["To"] = ", ".join(message.email_to)
        email_message["Subject"] = message.subject
        email_message.set_content(message.body, subtype="html")
        for attachment in message.attachments:
            attachment_maintype, attachment_subtype = attachment.content_type.split("/", 1)
            email_message.add_attachment(
                attachment.data,
                filename=attachment.filename,
                maintype=attachment_maintype,
                subtype=attachment_subtype
            )
        with smtplib.SMTP(self.smtp_settings.server, self.smtp_settings.port) as server:
            if self.smtp_settings.tls:
                server.starttls()
            if self.smtp_settings.username and self.smtp_settings.password:
                server.login(self.smtp_settings.username, self.smtp_settings.password)
            server.send_message(email_message)

    def _default_from_user(self):
        if self.smtp_settings.username:
            return self.smtp_settings.username
        return super()._default_from_user()

    def _default_from_domain(self):
        if self.smtp_settings.username and "@" in self.smtp_settings.username:
            return self.smtp_settings.username.split("@", 1)[1]
        return super()._default_from_domain()
