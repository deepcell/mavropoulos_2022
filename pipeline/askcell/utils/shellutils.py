r"""Infrastructure to assist in safely building and running complex command line applications.

:class:`Command` is a descendant of the builtin Python list class with a small amount of extra functionality that is
specialized toward storing command line arguments that are to be executed.  Its constructor and class instance has an
attribute named `shell` to differentiate which commands that are intended to be executed through the shell from those
that do not.  Its default value is set to False.  This behavioral distinction has implications for how a Command is
constructed, since list of arguments intended for shell execution require additional care.

This module is created to:

  1. unify the two disparate modes of command construction
  2. encourage use of non-shell execution
  3. make shell execution safer
  4. decouple the complexity of command creation from command execution.


Unify command construction
--------------------------

The Python subprocess module provides a variety of APIs that accept commands to execute either via:

  1. a list of arguments to be executed without a command shell by one of the subprocess APIs with shell=False.

  2. a string of properly escaped arguments concatenated with whitespace to be executed within a command shell by one of
     the subprocess APIs with shell=True.

:class:`Command` unifies these two styles of construction, by providing semantics that allow shell-bound commands to
remain as lists of arguments.

Encourage use of non-shell execution
------------------------------------

When possible, don't use shell execution.  Doing so makes your command impervious to errors that could lead to [Command
Injection]_ exploitations, since all arguments are processed literally and no shell-semantics are available.  e.g.  an
argument with an embedded space is treated as a single argument and not two arguments.

However, some commands benefit from functionality that a command shell provides and the extra complexity to ensure safe
execution is justified.  The examples of such cases include the command argument(s) with redirection, use of globbing
(wildcard characters), use of shell error handling, piping, and environment variables,

Make shell execution safer
--------------------------

If one must use shell execution, then all Command arguments are escaped using the :func:`shlex.quote` function except
those provided as :class:`StringLiteral` or :class:`ShellLiteral`.  So by default, all arguments are escaped and
treated as they would be if running with shell=False.

Decouple the complexity of command creation from command execution
------------------------------------------------------------------

Constructing command lines can be extremely complex.  A single command may have dozens of parameters, incorporate
complex shell semantics, and utilize multiple dedicated functions in its construction.  Likewise, command
execution can be complex.  Using the full :class:`subprocess.Popen` API often requires threading.

Module Description
==================

An executable command is represented by a single Command object, where :class:`ShellLiteral` and :class:`StringLiteral`,
a subclass of the builtin Python str class, is used to delineate arguments that are exempt from shell-escaping rules.

In this module, Command object with shell attribute set to True utilizes :class:`ShellLiteral` and calls
:func:`shell_string` to create a correctly constructed string argument for the shell.  The :class:`StringLiteral`
object is used to preserve the literal meaning of special characters.  An attempt to create a shell command string on
Command object with shell attribute set to False will throw an `NotImplementedError`.

Most of the use case for this module will be Command object in `shell=False` mode.  When `shell=False`, each string
element in Command object behaves as a :class:`Literal` object.

How to use this module:

   1. Build Command with `shell=False`, if possible.
   2. Build Command with `shell=True`
      (http://pubs.opengroup.org/onlinepubs/009695399/utilities/xcu_chap02.html#tag_02_07_04)

       * if command is expected to perform redirection, | & ; < >
       * if command contains shell commands such as echo, dirs, history
         (https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6/unix_tutorial/misc/unix_commands.html)
       * if command contains shell reserved words (!, do, esac, in, {, }, done, if, fi, elif, else, while, for,
         then, until, case)
       * if command is to handle the environment variable such as `export`
       * if command uses wildcard characters such as `*` and `?`

When to use :class:`StringLiteral` instances:

    1. Strings with characters that are pre-escaped, but do not require shell semantics.

When to use :class:`ShellLiteral` instances:

    1. Job control
    2. Variable replacement
    3. File redirection
    4. File globbing

Example of shell=False and its equivalence to shell=True::

    >>> cmd =  Command('sentieon', 'driver', '-i', 'in.bam', '--algo', 'Dedup', 'out.bam')              # Example 1
    >>> cmd2 = Command('sentieon', 'driver', '-i', 'in.bam', '--algo', 'Dedup', 'out.bam', shell=True)  # Example 2
    >>> cmd.command_list
    ['sentieon', 'driver', '-i', 'in.bam', '--algo', 'Dedup', 'out.bam']
    >>> cmd2
    ['sentieon', 'driver', '-i', 'in.bam', '--algo', 'Dedup', 'out.bam']
    >>> cmd2.shell_string
    'sentieon driver -i in.bam --algo Dedup out.bam'

Example for shell=True, run_cmd uses cmd.shell_string to construct command argument when shell=True::

    >>> cmd = Command(                                                                               # Example 3
    ...     'gzip', '-dc', 'in.vcf.gz',
    ...     ShellLiteral('|'),
    ...     'sentieon', 'util', 'vcfconvert', '-', 'in.vcf',
    ... )
    >>> cmd.shell_string
    'gzip -dc in.vcf.gz | sentieon util vcfconvert - in.vcf'

    >>> cmd = Command('ls', ShellLiteral('*'))                                                       # Example 4
    >>> cmd.shell_string
    'ls *'

    >>> cmd = Command('history', shell=True)                                                         # Example 5
    >>> cmd.shell_string
    'history'

    >>> cmd = Command('history', ShellLiteral('>'), '/mnt/scratch/temp.txt')                         # Example 6
    >>> cmd.shell_string
    'history > /mnt/scratch/temp.txt'

    >>> cmd = Command('export', 'SENTIEON_PYTHON=/bin/python', shell=True)                           # Example 7
    >>> cmd.shell_string
    'export SENTIEON_PYTHON=/bin/python'

    >>> cmd = Command('ls', 'anything')  # is equivalent to cmd = Command(['ls', 'anything'])        # Example 8
    >>> cmd
    ['ls', 'anything']


Example with intentionally silly semantics, when command contains wildcard character, yet Command object is not
correctly set as shell=True::

    >>> from .shellutils import run_cmd, Command

    >>> cmd = Command('ls', '*')                                                                     # Example 10
    >>> run_cmd(cmd)   # See Example 4 for correct use  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    CalledProcessError: Command '['ls', '*']' returned non-zero exit status 2.

Example with intentionally silly semantics, when shell-based command is used as Command object with shell=False::

    >>> cmd = Command('history')                                                                     # Example 11
    >>> run_cmd(cmd)  # See Example 5 for correct use
    Traceback (most recent call last):
    FileNotFoundError: [Errno 2] No such file or directory: 'history'

    >>> cmd = Command('export', 'SENTIEON_PYTHON=/bin/python')                                       # Example 12
    >>> run_cmd(cmd)  # See Example 7 for correct use
    Traceback (most recent call last):
    FileNotFoundError: [Errno 2] No such file or directory: 'export'


.. _[Command Injection]: https://www.kevinlondon.com/2015/07/26/dangerous-python-functions.html

"""

from __future__ import annotations

import shlex

from subprocess import check_call
from typing import Any, Iterable, Sequence, Union

from .pathutils import Path


CmdArg = Union[bool, int, str, float, "ShellLiteral", "StringLiteral", Path]
CmdList = list[CmdArg]


def split_cmd(cmd: str) -> list[str]:
    r"""Split strings using shell-like syntax.

    Args:
        cmd: command string

    Returns:
        list of command strings

    Example:
        >>> cmd = 'ls -l .'
        >>> split_cmd(cmd)
        ['ls', '-l', '.']

        >>> cmd = 'a b'
        >>> split_cmd(cmd)
        ['a', 'b']

        >>> cmd = "a b"
        >>> split_cmd(cmd)
        ['a', 'b']

        >>> cmd = 'a \\ b'
        >>> split_cmd(cmd)
        ['a', ' b']

    """
    return shlex.split(cmd)


def join_cmd(cmd: Sequence[str]) -> str:
    """Join command arguments with proper escaping and to protect from shell injection attacks.

    _Only_ a single command and its arguments should be passed into this function.
    The command and arguments are joined into a single string with proper escaping.
    This function is _not_ intended to take arbitrary shell commands,
    as the escaping will disable shell directives, such as redirection, job control,
    and other shell constructs.

    Good uses:

        >>> join_cmd(['ls', 'filename with space', '/path with semicolon;'])
        "ls 'filename with space' '/path with semicolon;'"

    Bad uses:

       # Don't run join_cmd on compound shell commands -- the shell controls will be escaped!
       >>> join_cmd(['ls', 'a b', '&&', 'true'])
       "ls 'a b' '&&' true"

    Instead:

        >>> cmd1 = join_cmd(['ls', 'a b'])
        >>> cmd2 = join_cmd(['true'])
        >>> cmd  = f'{cmd1} && {cmd2}'

    Args:
        cmd: list of command strings to join

    Returns:
        command string

    Example:
        >>> join_cmd(['ls', '-l', '.'])
        'ls -l .'

        >>> join_cmd(['a', 'b'])
        'a b'

        >>> join_cmd(['ls', '-l', 'file with space'])
        "ls -l 'file with space'"

        >>> join_cmd(['ls', '-l', 'file with > weirdcharacter'])
        "ls -l 'file with > weirdcharacter'"

        >>> join_cmd(['ls', '-al', '/mnt'])
        'ls -al /mnt'

        # Good use to construct compound command string
        >>> cmd1 = join_cmd(['ls', '-l', 'somefile'])
        >>> cmd2 = join_cmd(['rm', '-rf', '~'])
        >>> '; '.join([cmd1, cmd2])
        "ls -l somefile; rm -rf '~'"

        # Bad use, do not use `;` to construct compound command string
        >>> join_cmd(['ls', '-l', 'somefile; rm -rf ~'])
        "ls -l 'somefile; rm -rf ~'"

    """
    return " ".join(shlex.quote(arg) for arg in cmd)


class StringLiteral(str):
    """String literal class to indicate a string that should not be escaped when processed by :class:`Command`."""


class ShellLiteral(StringLiteral):
    """Shell literal clas to indicate a string that requires a shell when processed by :class:`Command`."""


def escape_args(cmd: Iterable[CmdArg]) -> list[str]:
    """Shell-escaped each argument except those that are :class:`Literal` instances, which are returned unmodified.

    Args:
        cmd: list for command arguments to be executed via a shell

    Returns:
        shell-escaped version of command arguments that can be safely used in shell-mode

    """
    return [arg if isinstance(arg, StringLiteral) else shlex.quote(str(arg)) for arg in cmd]


class Command(CmdList):
    """Command object of list type representing command for single execution.

    Args:
        cmd:   list of command elements to be pieced together to represent single execution
        shell: indicate whether object is to be executed by the shell or not. Default is set to False

    Attributes:
        shell: indicate whether object is to be executed by the shell or not. Default is set to False

    """

    __slots__ = ("shell",)

    shell: bool

    def __init__(self, *cmd: CmdArg | Iterable[CmdArg], shell: bool = False) -> None:
        """Initialize a Command."""
        # Command object expects list, however if multiple strings are passed in, it auto-converts into list object
        items = cmd[0] if len(cmd) == 1 and isinstance(cmd[0], list) else cmd

        super().__init__(items)

        self.shell = shell or any(isinstance(item, ShellLiteral) for item in items)

    def __eq__(self, other: object) -> bool:
        """Compare with another command.

        Args:
            other: another command

        Returns:
            True if equal to other, else False

        """
        if not isinstance(other, Command):
            return NotImplemented

        if self.shell != other.shell:
            return False

        if self.shell:
            return self.shell_string == other.shell_string
        else:
            return self.command_list == other.command_list

    @property
    def shell_string(self) -> str:
        """Construct a shell-mode safe string. This is used for Command obj whose shell mode is set to True."""
        if not self.shell:
            raise RuntimeError("This command is not designed to run via a shell for security reasons")

        return " ".join(escape_args(self))

    @property
    def command_list(self) -> list[str]:
        """Construct command list for command. This is used for Command obj whose shell mode is set to False."""
        if self.shell:
            raise RuntimeError("This command requires shell semantics and cannot be run via non-shell mode.")

        if any(isinstance(arg, ShellLiteral) for arg in self):
            raise RuntimeError("Shell literals are not valid to use in non-shell commands.")

        return [str(arg) for arg in self]

    def __add__(self, cmd: Iterable[CmdArg]) -> Command:
        """Combine Command objects.

        Args:
            cmd: External Command object to add to end of self.

        Returns:
            Extended Command object.

        """
        # Rename for clarity
        cmd1, cmd2 = self, cmd

        if not isinstance(cmd2, Command):
            cmd2 = Command(*cmd2)

        shell = cmd1.shell or cmd2.shell

        if shell:
            if not cmd1.shell:
                cmd1 = cmd1.convert_to_shell()
            if not cmd2.shell:
                cmd2 = cmd2.convert_to_shell()

        return Command(*cmd1, *cmd2, shell=shell)

    def __iadd__(self, cmd: Iterable[CmdArg]) -> Command:
        """Combine Command objects inplace.

        Args:
            cmd: External Command object to add to end of self.

        Returns:
            Extended Command object.

        """
        combined = self + cmd
        self[:] = combined
        self.shell = combined.shell
        return self

    def __and__(self, cmd: Command) -> Command:
        """Run self then cmd.

        Args:
            cmd: External Command object to add as secondary command to end of self.

        Returns:
            Extended Command object.

        """
        return self + Command(ShellLiteral(";")) + cmd

    def __or__(self, cmd: Command) -> Command:
        """Pipe output of self to cmd.

        Args:
            cmd: External Command object to pipe output of self to.

        Returns:
            Extended Command object.

        """
        return self + Command(ShellLiteral("|")) + cmd

    def redirect_output(
        self,
        *,
        stdout: CmdArg | None = None,
        stderr: CmdArg | None = None,
        append: bool = False,
    ) -> Command:
        """Redirect output streams.

        Args:
            stdout: standard output file
            stderr: standard error file
            append: True to append to output files, False to rewrite

        Returns:
            redirected command

        """
        cmd = self
        lt = ">>" if append else ">"

        if stdout is not None and stdout == stderr:
            cmd += Command(ShellLiteral(lt), stdout, ShellLiteral("2>&1"))
        else:
            if stdout is not None:
                cmd += Command(ShellLiteral(lt), stdout)

            if stderr is not None:
                cmd += Command(ShellLiteral(f"2{lt}"), stderr)

        return cmd

    def redirect_stdout_to_stderr(self) -> Command:
        """Redirect standard output to standard error.

        Returns:
            redirected command

        """
        return self + Command(ShellLiteral(">&2"))

    def redirect_stderr_to_stdout(self) -> Command:
        """Redirect standard error to standard outout.

        Returns:
            redirected command

        """
        return self + Command(ShellLiteral("2>&1"))

    def in_subshell(self) -> Command:
        """Run command in a subshell.

        Returns:
            subshell command

        """
        return Command(ShellLiteral("(")) + self + Command(ShellLiteral(")"))

    def convert_to_shell(self) -> Command:
        """Convert command to ensure Command shell attribute is True.

        A copy is returned if shell=True.

        Returns:
            new Command with shell=True.

        """
        return Command(self, shell=True)

    def as_output_pipe(self) -> Command:
        """Convert command to a named pipe using bash process substitution.

        Returns:
            new Command

        This function wraps a command to send its stdout to a named pipe.  For example:

        ...code: bash

            diff -u <(zcat file1.txt.gz) <(zcat file2.txt.gz)

        will allow two gzipped files to be compared using ``diff`` by uncompressing each file via a named pipe.

        See:
          * https://www.gnu.org/software/bash/manual/html_node/Process-Substitution.html
          * https://tldp.org/LDP/abs/html/process-sub.html

        """
        return Command(ShellLiteral("<(")) + self + Command(ShellLiteral(")"))


def run_cmd(
    cmd: str | list[str] | Command,
    strict_mode: bool = True,
    executable: str = "/bin/bash",
    **kwargs: Any,
) -> None:
    r"""Run a shell command with arguments using :py:func:`subprocess.check_call`.

    If ``cmd`` is not a Command, it is passed as is to :py:func:`subprocess.check_call`. By default the command
    will be run in non-shell mode unless ``shell=True`` is one of the keyworded arguments.

    If ``cmd`` is a Command and ``cmd.shell`` is False, it is run in non-shell mode.

    If ``cmd`` is a Command with ``cmd.shell`` set to True, it is run in shell mode. Given the default ``executable``,
    the shell will be changed from the Bourne shell (sh) to the Bourne Again shell (bash). If ``strict_mode`` ``cmd``
    will be run subject to `bash strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_ settings
    while keeping the default Internal Field Separator (IFS) (i.e. space, newline, and tab characters).

    Args:
        cmd: input command
        strict_mode: Allows running cmd in shell mode using "bash strict mode" settings
        executable: Changes shell mode executable from /bin/sh to /bin/bash (default)
        \*\*kwargs: dictionary of keyword arguments for :py:func:`subprocess.check_call`

    """
    if not isinstance(cmd, Command):
        # FIXME: Should this raise an error? or issue a warning?
        check_call(cmd, **kwargs)

    elif not cmd.shell:
        check_call(cmd.command_list, shell=False, **kwargs)

    else:
        shell_string = cmd.shell_string

        if strict_mode:
            if executable != "/bin/bash":
                raise ValueError(f"The shell executable must be /bin/bash if strict mode and not {executable}")

            shell_string = r"set -euo pipefail; IFS=$' \n\t'; " + shell_string

        check_call(shell_string, shell=True, executable=executable, **kwargs)
